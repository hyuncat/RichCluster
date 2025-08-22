#ifndef DavidClustering_h
#define DavidClustering_h

#include <Rcpp.h>
#include <string>
#include <vector>
#include <unordered_set>

class DavidClustering {
public:
    DavidClustering(
        const Rcpp::CharacterVector& terms,
        const Rcpp::CharacterVector& geneIDs,
        double similarityThreshold,
        int initialGroupMembership,
        int finalGroupMembership,
        double multipleLinkageThreshold
    );

    Rcpp::List run();

private:
    using TermSet = std::unordered_set<int>;

    void calculateKappaScores();
    void findInitialSeeds();
    void mergeSeeds();
    double calculateDiceCoefficient(const TermSet& seed1, const TermSet& seed2);

    // Input data
    Rcpp::CharacterVector terms;
    Rcpp::CharacterVector geneIDs;
    int n_terms;
    int totalGeneCount;

    // Parameters
    double similarityThreshold;
    int initialGroupMembership;
    int finalGroupMembership;
    double multipleLinkageThreshold;

    // Internal data structures
    std::vector<std::vector<double>> kappaMatrix;
    std::vector<TermSet> initialSeeds;
    std::vector<TermSet> finalClusters;
};

#endif // DavidClustering_h
