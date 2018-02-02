#include "ctane.h"
#include "../util/output.h"
#include "partitiontable.h"
#include <iterator>

CTane::CTane(Database& db)
    :BaseMiner(db) {
}

int CTane::nrCFDs() const {
    return fCFDs.count();
}

Itemset CTane::getConstantCFDs(const Itemset& nodeItems, const SimpleTidList& tids, const std::vector<MinerNode<PartitionTidList> >& attrs) {
    Itemset res;
    PartitionTidList nodeTidsPart = convert(tids);
    Itemset nodeAttrs = fDb.getAttrVector(nodeItems);
    for (const auto& attrNode : attrs) {
        int attr = -1 - attrNode.fItem;
        if (std::binary_search(nodeAttrs.begin(), nodeAttrs.end(), attr)) continue;
        auto isect = PartitionTable::intersection(nodeTidsPart, attrNode.fTids);
        if (isect.fNrSets == 1) {
            res.push_back(fDb.getRow(tids[0])[attr]);
        }
    }
    return res;
}

bool CTane::precedes(const Itemset& a, const Itemset& b) {
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

void CTane::pruneCands(std::vector<MinerNode<PartitionTidList> >& items, const Itemset& sub, int out) {
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

void CTane::mine(int minsup, double minconf, bool variable) {
    fMaxSize = fDb.nrAttrs();
    fMinSup = minsup;
    fMinConf = minconf;
    if (variable) {
        fAllAttrs = range(-fDb.nrAttrs(), 0);
    }
    auto allItems = range(1, fDb.nrItems()+1);
    fAllAttrs.insert(fAllAttrs.end(), allItems.begin(), allItems.end());
    PartitionTable::fDbSize = fDb.size();
    std::vector<MinerNode<PartitionTidList> > items;
    if (variable) {
        items = getPartitionSingletons();
    }
    else {
        items.clear();
    }
    fItemsLayer = getSingletons(fMinSup);
    auto constants = getItemsetLayer();
    for (const auto& cp : constants) {
        items.emplace_back(cp.first[0], convert(cp.second), cp.second.size());
    }
    for (auto& a : items) {
        Itemset at;
        for (int cat : fAllAttrs) {
            if (a.fItem < 0 && cat >= 0) continue;
            if (a.fItem != cat && fDb.getAttrIndex(a.fItem) == fDb.getAttrIndex(cat)) continue;
            at.push_back(cat);
        }
        a.fCands = at;
    }
    fCandStore = PrefixTree<Itemset, Itemset>();
    fGens.addMinGen(Itemset(), support(convert(iota(fDb.size()))), 0);
    fStore[Itemset()] = convert(iota(fDb.size()));
    fCandStore.insert(Itemset(), fAllAttrs);

    while (!items.empty()) {
        std::sort(items.begin(), items.end(), [](
                const MinerNode<PartitionTidList>& a, const MinerNode<PartitionTidList>& b) -> bool {
            return a.fPrefix < b.fPrefix || (a.fPrefix == b.fPrefix && a.fItem < b.fItem);
        });
        std::unordered_map<Itemset,PartitionTidList> newStore;
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
                auto storedSub = fStore.find(sub);
                if (storedSub == fStore.end()) {
                    continue;
                }
                if (out < 0) {
                    double e = PartitionTable::partitionError(storedSub->second, inode.fTids);
                    double conf = 1 - (e / support(storedSub->second));
                    if (conf >= fMinConf) {
                        fCFDs.addCFD(sub, out);
                    }
                    if (conf >= 1) {
                        inode.fCands = subset(inode.fCands, out);
                        inode.fCands = intersection(inode.fCands, iset);
                        pruneCands(items, sub, out);
                    }
                }
                else {
                    auto cleanVio = setdiff(storedSub->second.fTids, inode.fTids.fTids);
                    double conf = 1.0 - ((double)cleanVio.size() / support(storedSub->second));
                    if (conf >= fMinConf) {
                        fCFDs.addCFD(sub, out);
                    }
                    if (conf >= 1) {
                        inode.fCands = subset(inode.fCands, out);
                        inode.fCands = intersection(inode.fCands, iset);
                        pruneCands(items, sub, out);
                    }
                }
            }
            if (inode.fCands.empty()) continue;
            if (!fGens.addMinGen(where(iset, [](int i){return i >= 0;}), inode.fSupp, inode.fHash)) continue;

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
        fStore.swap(newStore);
        fCandStore = newCandStore;
        items.swap(suffix);
    }
}

std::vector<std::pair<Itemset,SimpleTidList> > CTane::getItemsetLayer() {
    std::vector<std::pair<Itemset,SimpleTidList> > res;
    std::vector<MinerNode<SimpleTidList> > suffix;
    for (int i = 0; i < fItemsLayer.size(); i++) {
        const auto &inode = fItemsLayer[i];
        Itemset iset = join(inode.fPrefix, inode.fItem);
        if (!fGens.addMinGen(iset, inode.fSupp, inode.fHash)) continue;

        res.emplace_back(iset, inode.fTids);
        for (int j = i + 1; j < fItemsLayer.size(); j++) {
            const auto &jnode = fItemsLayer[j];
            if (jnode.fPrefix != inode.fPrefix) continue;
            Itemset jset = join(jnode.fPrefix, jnode.fItem);
            Itemset newset = join(iset, jset);
            SimpleTidList ijtids = intersection(inode.fTids, jnode.fTids);
            int ijsupp = ijtids.size();
            if (ijsupp >= fMinSup && ijsupp != inode.fSupp){
                int jtem = newset.back();
                newset.pop_back();
                suffix.emplace_back(jtem, std::move(ijtids), ijsupp, newset);
            }
        }
    }
    std::sort(suffix.begin(), suffix.end());
    fItemsLayer.swap(suffix);
    return res;
}