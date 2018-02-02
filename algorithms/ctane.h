#ifndef ALGORITHMS_CTANE_H_
#define ALGORITHMS_CTANE_H_

#include "baseminer.h"

class CTane : public BaseMiner {
public:
    CTane(Database&);
    int nrCFDs() const;
    void mine(int, double=1, bool=true);
    void mineItemsets(const Itemset&, std::vector<MinerNode<SimpleTidList> >, const std::vector<MinerNode<PartitionTidList> >&);
    Itemset getConstantCFDs(const Itemset&, const SimpleTidList&, const std::vector<MinerNode<PartitionTidList> >&);
    void getCFDs(const Itemset&, const std::vector<MinerNode<PartitionTidList> >&);
    void getVariableCFDs(const Itemset&, std::vector<MinerNode<PartitionTidList> >&, std::vector<std::pair<Itemset, int> >&);
    void pruneCands(std::vector<MinerNode<PartitionTidList> >& items, const Itemset& sub, int out);
        std::vector<std::pair<Itemset,SimpleTidList> > getItemsetLayer();
    bool precedes(const Itemset& a, const Itemset& b);
private:
    GeneratorStore<int> fGens;
    std::unordered_map<Itemset,PartitionTidList> fStore;
    Itemset fAllAttrs;
    PrefixTree<Itemset, Itemset> fCandStore;
    std::vector<MinerNode<SimpleTidList> > fItemsLayer;
    double fMinConf;
};

#endif // ALGORITHMS_CTANE_H_