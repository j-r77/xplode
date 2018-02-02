#ifndef ALGORITHMS_XPLODE_H_
#define ALGORITHMS_XPLODE_H_

#include "baseminer.h"
#include <sstream>

struct RHSScoreNode {
    RHSScoreNode(int rhs, double score, double conf=0, double dirtyConf = 0, double bound = 0) {
        fRhs = rhs;
        fScore = score;
        fConf = conf;
        fBound = bound;
        fDirtyConf = dirtyConf;
    }

    bool operator< (const RHSScoreNode& rhs) const {
        return (fBound < rhs.fBound) || (fBound == rhs.fBound && fScore < rhs.fScore);
    }

    int fRhs;
    double fScore;
    double fConf;
    double fDirtyConf;
    double fBound;
};

struct FringeNode {
    FringeNode() {
        fExpanded = false;
    }
    FringeNode(const MinerNode<PartitionTidList>& node) {
        fNode = node;
        fExpanded = false;
        fRHSes = std::priority_queue<RHSScoreNode>();
    }
    bool operator< (const FringeNode& rhs) const {
        if (fRHSes.empty()) {
            if (rhs.fRHSes.empty()) {
                return fNode.fItem > rhs.fNode.fItem;
            }
            else {
                return false;
            }
        }
        else if (rhs.fRHSes.empty()) {
            return true;
        }
        return fRHSes.top() < rhs.fRHSes.top();
    }

    MinerNode<PartitionTidList> fNode;
    std::priority_queue<RHSScoreNode> fRHSes;
    bool fExpanded;
};

struct Modification {
    Modification(int t, int a, int d, int c) {
        fTid = t;
        fAttr = a;
        fVc = c;
        fVd = d;
    }
    int fTid;
    int fAttr;
    int fVc;
    int fVd;

    std::string toString() const {
        std::stringstream ss;
        ss << '{' << fTid << ',' << fAttr << ',' << fVc << ',' << fVd << '}';
        return ss.str();
    }

};

const unsigned ALL_SCORES = 1;
const unsigned FAST_TIDS = 2;
const unsigned CONVERT_EQ = 4;

class XPlode : public BaseMiner {
public:
    XPlode(Database&, Database&, unsigned optimize = ALL_SCORES | FAST_TIDS | CONVERT_EQ);
    XPlode(Database&, Database&, std::vector<int>);
    CFD explain(int, double=1, bool=true);
    CFD postExplain(int, double=1, bool=true);

    Itemset getDiffs();
    void computeModifications();
    void pruneCands(std::vector<MinerNode<PartitionTidList> >& items, const Itemset& sub, int out);
    std::vector<FringeNode> pruneFringe(const std::vector<FringeNode>& fringe, double score);
    bool precedes(const Itemset& a, const Itemset& b);
    std::vector<MinerNode<SimpleTidList> > getDirtySingletons(int minsup);
    std::vector<std::pair<Itemset,SimpleTidList> > getItemsetLayer();
    void mineDirtyItemsets();
    std::map<Itemset, int> convertCFD(const Itemset &lhs, int rhs, const PartitionTidList &tids);
    PartitionTidList pullback(const Itemset &iset, const PartitionTidList &cleanTids);
    PartitionTidList getTids(const Itemset &lhs, const Database &);
    std::vector<PartitionTidList> tidListIntersections(const PartitionTidList &lhs, const std::vector<PartitionTidList*> rhses);

    double score(const Itemset& lhs, int rhs, const SimpleTidList& vio);
    double upperBound(const Itemset &lhs, int rhs, const SimpleTidList &vio);

private:
    Database& fDirty;
    double fMinConf;
    unsigned fOptimize;

    std::map<int, PartitionTidList> fSingletons;
    GeneratorStore<int> fGeneratorCheck;
    SimpleTidList fRepairs;
    std::map<CFD, std::map<Itemset, int> > fCPMap;
    std::vector<std::vector<Modification> > fModifications;
    std::map<Itemset, SimpleTidList> fDirtySupports;
};

#endif // ALGORITHMS_XPLODE_H_