#include "xplode.h"
#include "../util/output.h"
#include "../data/databasereader.h"
#include <iterator>

XPlode::XPlode(Database& dirty, Database& rep, unsigned optimize)
        :BaseMiner(rep), fDirty(dirty), fOptimize(optimize) {
    fRepairs = Database::getDiffs(dirty, rep);
}

XPlode::XPlode(Database& dirty, Database& rep, std::vector<int> repairs)
        :BaseMiner(rep), fDirty(dirty) {
    fRepairs = repairs;
}

Itemset XPlode::getDiffs() {
    // Create set diffs of all items A and B
    // that occur in a modification A->B
    Itemset diffs;
    std::vector<int> attrMap(fDb.nrAttrs());
    for (const auto& mods : fModifications) {
        for (const auto& m : mods) {
            attrMap[m.fAttr] = 1;
        }
    }
    for (int ix = 0; ix < attrMap.size(); ix++) {
        if (attrMap[ix]) diffs.push_back(ix);
    }
    return diffs;
}

bool XPlode::precedes(const Itemset& a, const Itemset& b) {
    if (a.size() > b.size()) return false;
    if (a == b) return false;
    auto bAttrs = fDb.getAttrVector(b);
    for (int i : a) {
        if (!std::binary_search(b.begin(), b.end(), i)) {
            if (i < 0 && !std::binary_search(bAttrs.begin(), bAttrs.end(), -1-i)) {
                return false;
            }
            else if (i >= 0) {
                return false;
            }
        }
    }
    return true;
}

void XPlode::pruneCands(std::vector<MinerNode<PartitionTidList> >& items, const Itemset& sub, int out) {
    Itemset cset = join(sub, out);
    for (int i = 0; i < items.size(); i++) {
        MinerNode<PartitionTidList> &inode = items[i];
        const Itemset iset = join(inode.fPrefix, inode.fItem);
        if (precedes(sub, iset)) {
            if (std::find(inode.fCands.begin(), inode.fCands.end(), out) != inode.fCands.end()) {
                inode.fCands = subset(inode.fCands, out);
            }
            if (precedes(cset, iset)) {
                inode.fCands = intersection(inode.fCands, sub);
            }
        }
    }
}

PartitionTidList XPlode::getTids(const Itemset &iset, const Database &) {
    if (!(fOptimize & FAST_TIDS)) {
        return getTids(iset, fDb);
    }
    PartitionTidList res { SimpleTidList(), 0 };
    Itemset tp = where(iset, [](int i){return i >= 0;});
    Itemset vp = where(iset, [](int i){return i < 0;});
    bool first = true;
    for (int c : tp) {
        if (first) {
            res = fSingletons[c];
            first = false;
        }
        else {
            res = convert(intersection(res.fTids, fSingletons[c].fTids));
        }
        if (!support(res)) return res;
    }
    for (int v : vp) {
        if (first) {
            res = fSingletons[v];
            first = false;
        }
        else {
            res = PartitionTable::intersection(res, fSingletons[v]);
        }
    }
    return res;
}

std::vector<FringeNode> XPlode::pruneFringe(const std::vector<FringeNode>& fringe, double score) {
    std::vector<FringeNode> newFringe;
    for (auto &tmp : fringe) {
        if (tmp.fRHSes.size() && tmp.fRHSes.top().fBound > score)
            newFringe.push_back(tmp);
    }
    std::make_heap(newFringe.begin(), newFringe.end());
    return newFringe;
}

CFD XPlode::explain(int minsup, double minconf, bool variable) {
    fDb.toFront(fRepairs);
    fDirty.toFront(fRepairs);
    fRepairs = iota(fRepairs.size());
    fMinSup = minsup;
    fMinConf = minconf;
    PartitionTable::fDbSize = fDb.size();
    std::map<Itemset,PartitionTidList> tidmap;
    std::map<Itemset,PartitionTidList> dirtyTidmap;
    PrefixTree<Itemset, Itemset> candidateMap;
    Itemset allItems;

    computeModifications();
    mineDirtyItemsets();

    std::vector<MinerNode<PartitionTidList> > items;
    if (variable) {
        allItems = range(-fDb.nrAttrs(), 0);
        items = getPartitionSingletons();
    }

    auto constants = getSingletons(fMinSup);
    for (const auto& cp : constants) {
        allItems.push_back(cp.fItem);
        items.emplace_back(cp.fItem, convert(cp.fTids), cp.fTids.size());
    }
    std::sort(allItems.begin(), allItems.end());

    candidateMap = PrefixTree<Itemset, Itemset>();
    fGeneratorCheck.addMinGen(Itemset(), support(convert(iota(fDb.size()))), 0);
    tidmap[Itemset()] = convert(iota(fDb.size()));
    candidateMap.insert(Itemset(), allItems);
    Itemset diffs = getDiffs();
    std::vector<FringeNode> fringe;

    for (auto& a : items) {
        Itemset at;
        for (int cat : allItems) {
            if (a.fItem < 0 && cat >= 0) continue;
            if (a.fItem != cat && fDb.getAttrIndex(a.fItem) == fDb.getAttrIndex(cat)) continue;
            at.push_back(cat);
        }
        a.fCands = at;
        double bound = a.fItem < 0 ? fRepairs.size() : intersection(a.fTids.fTids, fRepairs).size();
        bound -= ((double)join(diffs, itemset(fDb.getAttrIndex(a.fItem))).size()/1000.0);
        if (bound > 0) {
            FringeNode fn(a);
            fn.fRHSes.emplace(a.fItem, 0, 0, 0, bound);
            tidmap[itemset(a.fItem)] = a.fTids;
            candidateMap.insert(itemset(a.fItem), a.fCands);
            fSingletons[a.fItem] = a.fTids;
            fringe.push_back(fn);
        }
    }
    std::make_heap(fringe.begin(), fringe.end());
    std::set<Itemset> visited;
    double maxScore = 0;
    double maxFoundScore = 0;
    CFD maxCfd;

    while (!fringe.empty()) {
        std::pop_heap(fringe.begin(), fringe.end());
        FringeNode top = fringe.back();
        fringe.pop_back();
        const Itemset iset = join(top.fNode.fPrefix, top.fNode.fItem);
        auto nodeAttrs = fDb.getAttrVector(iset);
        if (top.fRHSes.size() && top.fRHSes.top().fScore) {

            auto rhsNode = top.fRHSes.top();
            if (rhsNode.fBound <= maxFoundScore) break;
            top.fRHSes.pop();
            Itemset sub = subset(iset, rhsNode.fRhs);
            int out = rhsNode.fRhs;

            if (std::binary_search(top.fNode.fCands.begin(), top.fNode.fCands.end(), out)) {
                if (rhsNode.fConf >= fMinConf && rhsNode.fDirtyConf < rhsNode.fConf) {
                    if (rhsNode.fScore > maxFoundScore) {
                        maxFoundScore = rhsNode.fScore;
                        maxCfd = std::make_pair(sub, out);
                        if (rhsNode.fScore == rhsNode.fBound) break;
                    }
                }
            }
        }
        if (top.fRHSes.size() && top.fRHSes.top().fScore && top.fRHSes.top().fBound > maxScore) {
            fringe.push_back(top);
            std::push_heap(fringe.begin(), fringe.end());
        }

        if (top.fExpanded) continue;
        top.fExpanded = true;

        if (!fGeneratorCheck.addMinGen(where(iset, [](int i){return i >= 0;}), top.fNode.fSupp, top.fNode.fHash)) {
            top.fNode.fCands.clear();
        }

        if (top.fNode.fCands.size()) {
            tidmap[iset] = top.fNode.fTids;
        }
        else {
            tidmap[iset] = {SimpleTidList(), 0};
            continue;
        }


        auto extensions = setdiff(allItems, iset);
        std::vector<PartitionTidList*> expands;
        std::vector<MinerNode<PartitionTidList> > tmpSuffix;
        for (int e : extensions) {
            if (std::binary_search(nodeAttrs.begin(), nodeAttrs.end(), fDb.getAttrIndex(e))) continue;
            auto rsub = join(top.fNode.fPrefix, e);

            Itemset newset = join(iset, e);
            auto insSucc = visited.insert(newset);
            if (!insSucc.second) continue;

            auto c = top.fNode.fCands;

            for (int zz : newset) {
                const Itemset* subCands = candidateMap.find(subset(newset, zz));
                if (subCands) {
                    c = intersection(c, *subCands);
                }
                if (!c.size()) break;
            }
            candidateMap.insert(newset, c);
            if (c.size()) {
                auto subTidsPtr = tidmap.find(rsub);
                PartitionTidList subTids;
                if (subTidsPtr == tidmap.end()) {
                    auto insPtr = tidmap.emplace(rsub, getTids(rsub, fDb));
                    subTidsPtr = insPtr.first;
                }
                expands.push_back(&subTidsPtr->second);
                tmpSuffix.emplace_back(newset.back());
                tmpSuffix.back().fCands = c;
                newset.pop_back();
                tmpSuffix.back().fPrefix = newset;
            }
        }

        std::vector<MinerNode<PartitionTidList> > suffix;
        std::vector<PartitionTidList> exps;
        if (variable)
            exps = PartitionTable::intersection(top.fNode.fTids, expands);
        else
            exps = tidListIntersections(top.fNode.fTids, expands);
        for (int e = 0; e < exps.size(); e++) {
            if (support(exps[e]) >= fMinSup) {
                suffix.emplace_back(tmpSuffix[e].fItem, exps[e]);
                suffix.back().fCands = tmpSuffix[e].fCands;
                suffix.back().fPrefix = tmpSuffix[e].fPrefix;
                FringeNode fn(suffix.back());
                const Itemset jset = join(suffix.back().fPrefix, suffix.back().fItem);

                auto insect = intersection(jset, suffix.back().fCands);
                auto rhses = Itemset();
                for (int out : insect) {
                    Itemset sub = subset(jset, out);
                    if (!isValid(sub, out)) {
                        continue;
                    }

                    auto subTidsPtr = tidmap.find(sub);
                    PartitionTidList subTids;
                    if (subTidsPtr != tidmap.end()) {
                        subTids = subTidsPtr->second;
                    }
                    else {
                        subTids = getTids(sub, fDb);
                        tidmap[sub] = subTids;
                    }
                    if (support(subTids) <= 0) continue;
                    auto dSubTidsPtr = dirtyTidmap.find(sub);
                    PartitionTidList dSubTids;
                    if (dSubTidsPtr != dirtyTidmap.end()) {
                        dSubTids = dSubTidsPtr->second;
                    }
                    else {
                        dSubTids = pullback(sub, subTids);
                        dirtyTidmap[sub] = dSubTids;
                    }
                    auto dITidsPtr = dirtyTidmap.find(jset);
                    PartitionTidList dITids;
                    if (dITidsPtr != dirtyTidmap.end()) {
                        dITids = dITidsPtr->second;
                    }
                    else {
                        dITids = pullback(jset, suffix.back().fTids);
                        dirtyTidmap[jset] = dITids;
                    }
                    auto vios = out < 0 ? PartitionTable::violations(dSubTids, dITids) : setdiff(dSubTids.fTids, dITids.fTids);
                    double bound = upperBound(sub, out, vios);
                    if (bound > maxScore) {
                        double conf = 0;
                        double scor = 0;
                        double dirtyConf = 0;
                        bool notVioInCleaned;
                        if (out < 0) {
                            notVioInCleaned = !PartitionTable::violatedInCleaned(subTids, suffix.back().fTids, fRepairs.size());
                        }
                        else {
                            auto cleanVio = setdiff(subTids.fTids, suffix.back().fTids.fTids);
                            notVioInCleaned = (cleanVio.empty() || cleanVio[0] > fRepairs.back());
                        }
                        if (notVioInCleaned) {
                            double e = out < 0 ? PartitionTable::partitionError(subTids, suffix.back().fTids) : setdiff(subTids.fTids, suffix.back().fTids.fTids).size();
                            conf = 1.0 - (e / (double) support(subTids));
                            if (conf >= fMinConf) {
                                if (out < 0) {
                                    double dirtyE = PartitionTable::partitionError(dSubTids, dITids);
                                    dirtyConf = 1 - (dirtyE / support(dSubTids));
                                }
                                else {
                                    dirtyConf = 1.0 - (vios.size() / (double) support(subTids));
                                }
                                if (!(fOptimize & ALL_SCORES)) {
                                    if (conf < dirtyConf) {
                                        if (out < 0) {
                                            fCPMap[std::make_pair(sub, out)] = convertCFD(sub, out, suffix.back().fTids);
                                        }
                                        else {
                                            auto sSub = sub;
                                            std::sort(sSub.begin(), sSub.end(),
                                                      [this](const int & a, const int & b) -> bool
                                                      {
                                                          return this->fDb.getAttrIndex(a) < this->fDb.getAttrIndex(b) ;
                                                      });
                                            fCPMap[std::make_pair(sub, out)].emplace(sSub, out);
                                        }
                                        scor = score(sub, out, vios);
                                        if (scor > maxScore && conf >= fMinConf && dirtyConf < conf) {
                                            maxScore = scor;
                                            fringe = pruneFringe(fringe, maxScore);
                                        }
                                    }
                                }
                                if (conf >= 1) {
                                    top.fNode.fCands = intersection(top.fNode.fCands, subset(iset,out));
                                }
                            }
                        }
                        if (fOptimize & ALL_SCORES) {
                            if (out < 0) {
                                fCPMap[std::make_pair(sub, out)] = convertCFD(sub, out, suffix.back().fTids);
                            } else {
                                auto sSub = sub;
                                std::sort(sSub.begin(), sSub.end(),
                                          [this](const int &a, const int &b) -> bool {
                                              return this->fDb.getAttrIndex(a) < this->fDb.getAttrIndex(b);
                                          });
                                fCPMap[std::make_pair(sub, out)].emplace(sSub, out);
                            }
                            scor = score(sub, out, vios);
                            if (scor > maxScore && conf >= fMinConf && dirtyConf < conf) {
                                maxScore = scor;
                                fringe = pruneFringe(fringe, maxScore);
                            }
                        }
                        fn.fRHSes.emplace(out, scor, conf, dirtyConf, bound);
                    }
                    else {
                        suffix.back().fCands = subset(suffix.back().fCands, out);
                    }
                }
                fn.fNode.fCands = join(suffix.back().fCands, rhses);
                if (fn.fRHSes.size() && fn.fNode.fCands.size()) {
                    fringe.push_back(fn);
                }
                candidateMap.insert(jset, suffix.back().fCands);
            }
        }
        std::make_heap(fringe.begin(), fringe.end());
    }
    return maxCfd;
}

PartitionTidList XPlode::pullback(const Itemset &iset, const PartitionTidList &cleanTids) {
    if (!(fOptimize & CONVERT_EQ)) {
        return getTids(iset, fDirty);
    }
    std::map<Itemset, Itemset> lookup;
    auto nodeAttrs = fDb.getAttrVector(iset);
    auto tp = where(iset, [](int i){return i >= 0;});
    if (tp.size() == iset.size()) {
        return PartitionTidList(join(fDirtySupports[tp], setdiff(cleanTids.fTids, fRepairs)), 1);
    }
    SimpleTidList toMove = fDirtySupports[tp];
    for (int i : toMove) {
        auto vals = projection(fDirty.getRow(i), nodeAttrs);
        lookup[vals].push_back(i);
    }
    PartitionTidList res;
    res.fNrSets = 1;
    res.fTids.reserve(cleanTids.fTids.size());
    bool first = true;
    bool inserted = false;
    for (int ct : cleanTids.fTids) {
        if (ct == PartitionTidList::SEP) {
            if (inserted) {
                res.fTids.push_back(ct);
            }
            res.fNrSets++;
            first = true;
            inserted = false;
        }
        else if (ct > fRepairs.back()) {
            inserted = true;
            if (first) {
                first = false;
                auto vals = projection(fDirty.getRow(ct), nodeAttrs);
                if (lookup.find(vals) != lookup.end()) {
                    auto &ins = lookup[vals];
                    res.fTids.insert(res.fTids.end(), ins.begin(), ins.end());
                    ins.clear();
                }
            }
            res.fTids.push_back(ct);
        }
    }
    for (const auto& kvp : lookup) {
        if (kvp.second.size()) {
            if (res.fTids.size()) {
                res.fTids.push_back(PartitionTidList::SEP);
                res.fNrSets++;
            }
            res.fTids.insert(res.fTids.end(), kvp.second.begin(), kvp.second.end());
        }
    }
    return res;
}

std::vector<PartitionTidList> XPlode::tidListIntersections(const PartitionTidList &lhs, const std::vector<PartitionTidList*> rhses) {
    std::vector<PartitionTidList> res;
    res.reserve(rhses.size());
    for (const PartitionTidList* rhs : rhses) {
        res.emplace_back(intersection(lhs.fTids, rhs->fTids), 0 );
    }
    return res;
}

CFD XPlode::postExplain(int minsup, double minconf, bool variable) {
    fDb.toFront(fRepairs);
    fDirty.toFront(fRepairs);

    fRepairs = iota(fRepairs.size());
    fMaxSize = fDb.nrAttrs();
    fMinSup = minsup;
    fMinConf = minconf;
    CFDList globalExplanations;
    std::map<Itemset,PartitionTidList> tidmap;
    std::map<Itemset,PartitionTidList> dirtyTidmap;
    std::unordered_map<CFD, SimpleTidList> violations;
    PrefixTree<Itemset, Itemset> candidateMap;

    Itemset allItems;
    if (variable) {
        allItems = range(-fDb.nrAttrs(), 0);
    }
    PartitionTable::fDbSize = fDb.size();
    std::vector<MinerNode<PartitionTidList> > items;
    if (variable) {
        items = getPartitionSingletons();
    }
    else {
        items.clear();
    }
    mineDirtyItemsets();
    auto constants = getSingletons(fMinSup);
    for (const auto& cp : constants) {
        allItems.push_back(cp.fItem);
        items.emplace_back(cp.fItem, convert(cp.fTids), cp.fTids.size());
    }
    std::sort(allItems.begin(), allItems.end());
    for (auto& a : items) {
        Itemset at;
        for (int cat : allItems) {
            if (a.fItem < 0 && cat >= 0) continue;
            if (a.fItem != cat && fDb.getAttrIndex(a.fItem) == fDb.getAttrIndex(cat)) continue;
            at.push_back(cat);
        }
        a.fCands = at;
    }
    candidateMap = PrefixTree<Itemset, Itemset>();
    fGeneratorCheck.addMinGen(Itemset(), support(convert(iota(fDb.size()))), 0);
    tidmap[Itemset()] = convert(iota(fDb.size()));
    candidateMap.insert(Itemset(), allItems);
    while (!items.empty()) {
        std::sort(items.begin(), items.end(), [](
                const MinerNode<PartitionTidList>& a, const MinerNode<PartitionTidList>& b) -> bool {
            return a.fPrefix < b.fPrefix || (a.fPrefix == b.fPrefix && a.fItem < b.fItem);
        });
        std::map<Itemset,PartitionTidList> newStore;
        PrefixTree<Itemset, Itemset> newCandStore;
        std::vector<std::vector<std::pair<Itemset,int> > > next(items.size());
        for (int i = 0; i < items.size(); i++) {
            MinerNode<PartitionTidList>& inode = items[i];
            const Itemset iset = join(inode.fPrefix, inode.fItem);
            auto insect = intersection(iset, inode.fCands);

            for (int out : insect) {
                Itemset sub = subset(iset, out);
                // Discard Variable -> Constant and Constant -> Variable CFDs
                if (out < 0) {
                    if (sub.size() && !has(sub, [](int si) -> bool { return si < 0; })) continue;
                }
                else {
                    if (!sub.size() || has(sub, [](int si) -> bool { return si < 0; })) continue;
                }
                auto storedSub = tidmap.find(sub);
                if (storedSub == tidmap.end()) {
                    continue;
                }
                double conf;
                if (out < 0) {
                    double e = PartitionTable::partitionError(storedSub->second, inode.fTids);
                    conf = 1 - (e / support(storedSub->second));
                    if (conf >= fMinConf) {
                        if (!PartitionTable::violatedInCleaned(storedSub->second, inode.fTids, fRepairs.size())) {
                            const auto cfd = std::make_pair(sub, out);
                            const auto dirtySub = pullback(sub, storedSub->second);
                            const auto dirtySet = pullback(iset, inode.fTids);
                            auto vios = PartitionTable::violations(dirtySub, dirtySet);
                            if (vios.size() && vios[0] <= fRepairs.back()) {
                                double dirtyE = PartitionTable::partitionError(dirtySub, dirtySet);
                                double dirtyConf = 1 - (dirtyE / support(dirtySub));
                                if (dirtyConf < conf-(1.0/fDb.size())) {
                                    globalExplanations.emplace_back(sub, out);
                                    violations[cfd] = vios;
                                    fCPMap[cfd] = convertCFD(sub, out, inode.fTids);
                                }
                            }
                        }
                    }
                }
                else {
                    auto cleanVio = setdiff(storedSub->second.fTids, inode.fTids.fTids);
                    conf = 1.0 - ((double)cleanVio.size() / support(storedSub->second));
                    if (conf >= fMinConf) {
                        if (cleanVio.empty() || cleanVio[0] > fRepairs.back()) {
                            const auto cfd = std::make_pair(sub, out);
                            const auto dirtySub = pullback(sub, storedSub->second);
                            const auto dirtySet = pullback(iset, inode.fTids);
                            auto vios = setdiff(dirtySub.fTids, dirtySet.fTids);
                            if (vios.size() && vios[0] <= fRepairs.back()) {
                                double dirtyConf = 1 - ((double) vios.size() / (double) support(dirtySub));
                                if (dirtyConf < conf) {
                                    globalExplanations.emplace_back(sub, out);
                                    violations[cfd] = vios;
                                    auto sSub = sub;
                                    std::sort(sSub.begin(), sSub.end(),
                                              [this](const int & a, const int & b) -> bool
                                              {
                                                  return this->fDb.getAttrIndex(a) < this->fDb.getAttrIndex(b) ;
                                              });
                                    fCPMap[std::make_pair(sub, out)].emplace(sSub, out);
                                }
                            }
                        }
                    }
                }
                if (conf >= 1) {
                    inode.fCands = intersection(inode.fCands, subset(iset, out));
                    pruneCands(items, sub, out);
                }
            }
            if (inode.fCands.empty()) continue;
            if (!fGeneratorCheck.addMinGen(where(iset, [](int i){return i >= 0;}), inode.fSupp, inode.fHash)) continue;

            newStore[iset] = inode.fTids;
            auto nodeAttrs = fDb.getAttrVector(iset);
            for (int j = i+1; j < items.size(); j++) {
                if (j == i) continue;
                const auto& jnode = items[j];
                if (jnode.fPrefix != inode.fPrefix) continue;
                if (std::binary_search(nodeAttrs.begin(), nodeAttrs.end(), fDb.getAttrIndex(jnode.fItem))) continue;
                next[i].emplace_back(iset, j);
            }
        }
        for (int i = 0; i < items.size(); i++) {
            MinerNode<PartitionTidList> &inode = items[i];
            const Itemset iset = join(inode.fPrefix, inode.fItem);
            newCandStore.insert(iset, inode.fCands);
        }

        std::vector<MinerNode<PartitionTidList> > suffix;
        for (int i = 0; i < items.size(); i++) {
            std::vector<PartitionTidList*> expands;
            std::vector<MinerNode<PartitionTidList> > tmpSuffix;
            for (auto& newsetTup : next[i]) {
                int j = newsetTup.second;
                Itemset newset = join(newsetTup.first, items[j].fItem);
                auto c = items[i].fCands;
                for (int zz : newset) {
                    auto zsub = subset(newset, zz);
                    auto storedSub = newStore.find(zsub);
                    if (storedSub == newStore.end()) {
                        c.clear();
                        break;
                    }
                    const Itemset &subCands = *newCandStore.find(zsub);
                    c = intersection(c, subCands);
                }
                if (c.size()) {
                    expands.push_back(&items[j].fTids);
                    tmpSuffix.emplace_back(items[j].fItem);
                    tmpSuffix.back().fCands = c;
                    tmpSuffix.back().fPrefix = newsetTup.first;
                }
            }
            const auto exps = PartitionTable::intersection(items[i].fTids, expands);
            for (int e = 0; e < exps.size(); e++) {
                if (support(exps[e]) >= fMinSup) {
                    suffix.emplace_back(tmpSuffix[e].fItem, exps[e]);
                    suffix.back().fCands = tmpSuffix[e].fCands;
                    suffix.back().fPrefix = tmpSuffix[e].fPrefix;
                }
            }
        }
        tidmap.swap(newStore);
        candidateMap = newCandStore;
        items.swap(suffix);
    }

    computeModifications();

    std::sort(globalExplanations.begin(), globalExplanations.end(),
              [this,violations](const CFD &a, const CFD &b) -> bool {
                  auto ascore = score(a.first, a.second, violations.at(a));
                  auto bscore = score(b.first, b.second, violations.at(b));
                  return (ascore < bscore) || (ascore == bscore && a.second < b.second);
              });
    return globalExplanations.back();
}

std::map<Itemset, int> XPlode::convertCFD(const Itemset &lhs, int rhs, const PartitionTidList &tids) {
    std::map<Itemset, int> res;
    auto lhsAttrs = fDb.getAttrVector(lhs);
    int rAttr = -1-rhs;
    std::map<Itemset, std::map<int,int> > rules;
    int count = 0;
    for (int pi = 0; pi <= tids.fTids.size(); pi++) {
        if (pi == tids.fTids.size() || tids.fTids[pi] == PartitionTidList::SEP ) {
            const Transaction& trans = fDb.getRow(tids.fTids[pi-1]);
            Itemset lhsConstants = projection(trans, lhsAttrs);
            rules[lhsConstants][trans[rAttr]] = count;
            count = 0;
        }
        else {
            count++;
        }
    }
    for (const auto& kvp : rules) {
        int max = 0;
        CFD maxCFD;
        std::map<Itemset, int> cfds;
        for (const auto& rhses : kvp.second) {
            if (rhses.second > max) {
                max = rhses.second;
                maxCFD = std::make_pair(kvp.first, rhses.first);
            }
        }
        res[maxCFD.first] = maxCFD.second;
    }
    return res;
}

double XPlode::score(const Itemset& lhs, int rhs, const SimpleTidList& vios) {
    double c = 0;
    auto lhsAttrs = fDb.getAttrVector(lhs);
    const auto& cps = fCPMap.at(std::make_pair(lhs, rhs));
    std::vector<int> tp(fDb.nrAttrs());
    for (int item : lhs) {
        tp[fDb.getAttrIndex(item)] = item;
    }
    int rhsA = fDb.getAttrIndex(rhs);
    tp[rhsA] = rhs;
    for (int t : vios) {
        if (t <= fRepairs.back()) {
            for (const auto& m : fModifications[t]) {
                auto lhsValues2 = projection(fDirty.getRow(m.fTid), lhsAttrs);
                if (m.fAttr == rhsA) {
                    auto lhsValues = projection(fDirty.getRow(m.fTid), lhsAttrs);
                    if (cps.find(lhsValues) != cps.end()) {
                        if (m.fVc == cps.at(lhsValues)) {
                            c++;
                        }
                    }
                } else if (std::find(lhsAttrs.begin(), lhsAttrs.end(), m.fAttr) != lhsAttrs.end()) {
                    if (tp[m.fAttr] >= 0) {
                        c++;
                    } else {
                        auto lhsValues = projection(fDb.getRow(m.fTid), lhsAttrs);
                        if (cps.find(lhsValues) != cps.end()) {
                            if (fDb.getRow(m.fTid)[rhsA] == cps.at(lhsValues)) {
                                c++;
                            }
                        }
                    }
                }
            }
        }
        else {
            break;
        }
    }
    return c - ((lhs.size()+1.0)/(double)fDb.nrAttrs()+1.0);
}

double XPlode::upperBound(const Itemset &lhs, int rhs, const SimpleTidList &vios) {
    double c = 0;
    auto attrs = fDb.getAttrVector(join(lhs, rhs));
    std::vector<int> modAttrMask(fDb.nrAttrs());
    for (int t : vios) {
        if (t <= fRepairs.back()) {
            for (const auto &m : fModifications[t]) {
                modAttrMask[m.fAttr] = 1;
                c++;
            }
        } else {
            break;
        }
    }
    std::vector<int> modAttrs;
    for (int i = 0; i < modAttrMask.size(); i++) {
        modAttrs.push_back(i);
    }
    std::vector<int> neededAttrs = join(modAttrs, attrs);
    return c - (neededAttrs.size() / 1000.0);
}

std::vector<MinerNode<SimpleTidList> > XPlode::getDirtySingletons(int minsup) {
    std::vector<MinerNode<SimpleTidList> > singletons;
    std::unordered_map<int, SimpleTidList> nodeIndices;
    for (unsigned item = 1; item <= fDb.nrItems(); item++) {
        if (fDb.frequency(item) >= minsup) {
            nodeIndices[item] = SimpleTidList();
        }
    }

    for (unsigned row : fRepairs) {
        const auto& tup = fDirty.getRow(row);
        for (int item : tup) {
            const auto& it = nodeIndices.find(item);
            if (it != nodeIndices.end()) {
                it->second.push_back(row);
            }
        }
    }

    for (const auto& kvp : nodeIndices) {
        singletons.push_back(MinerNode<SimpleTidList>(kvp.first, kvp.second));
        singletons.back().hashTids();
    }

    std::sort(singletons.begin(), singletons.end());
    return singletons;
}

void XPlode::mineDirtyItemsets() {
    std::vector<MinerNode<SimpleTidList> > items = getDirtySingletons(fMinSup);
    fDirtySupports[Itemset()] = fRepairs;
    while (!items.empty()) {
        std::vector<MinerNode<SimpleTidList> > suffix;
        for (int i = 0; i < items.size(); i++) {
            const auto &inode = items[i];
            Itemset iset = join(inode.fPrefix, inode.fItem);
            fDirtySupports[iset] = inode.fTids;
            for (int j = i + 1; j < items.size(); j++) {
                const auto &jnode = items[j];
                if (jnode.fPrefix != inode.fPrefix) continue;
                Itemset jset = join(jnode.fPrefix, jnode.fItem);
                Itemset newset = join(iset, jset);
                SimpleTidList ijtids = intersection(inode.fTids, jnode.fTids);
                int ijsupp = ijtids.size();
                if (ijsupp) {
                    int jtem = newset.back();
                    newset.pop_back();
                    suffix.emplace_back(jtem, std::move(ijtids), ijsupp, newset);
                }
            }
        }
        std::sort(suffix.begin(), suffix.end());
        items.swap(suffix);
    }
}

void XPlode::computeModifications() {
    fModifications.clear();
    fModifications.resize(fRepairs.size());
    for (int i : fRepairs) {
        fModifications[i] = std::vector<Modification>();
        const auto &tup = fDirty.getRow(i);
        const auto &rep = fDb.getRow(i);
        for (int a = 0; a < tup.size(); a++) {
            if (tup[a] != rep[a]) {
                fModifications[i].emplace_back(i, a, tup[a], rep[a]);
            }
        }
    }
}