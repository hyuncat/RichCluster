#include "DavidClustering.h"
#include "StringUtils.h"
#include <Rcpp.h>

DavidClustering::DavidClustering(
    const Rcpp::CharacterVector& terms,
    const Rcpp::CharacterVector& geneIDs,
    double similarityThreshold,
    int initialGroupMembership,
    int finalGroupMembership,
    double multipleLinkageThreshold
) : terms(terms),
    geneIDs(geneIDs),
    similarityThreshold(similarityThreshold),
    initialGroupMembership(initialGroupMembership),
    finalGroupMembership(finalGroupMembership),
    multipleLinkageThreshold(multipleLinkageThreshold) {

    n_terms = terms.size();
    totalGeneCount = StringUtils::countUniqueElements(geneIDs);
    kappaMatrix.resize(n_terms, std::vector<double>(n_terms, 0.0));
}

Rcpp::List DavidClustering::run() {
    Rcpp::Rcout << "Calculating kappa scores..." << std::endl;
    calculateKappaScores();

    Rcpp::Rcout << "Finding initial seeds..." << std::endl;
    findInitialSeeds();

    Rcpp::Rcout << "Merging seeds..." << std::endl;
    mergeSeeds();

    // Format final clusters for output
    std::vector<int> clusterColumn;
    std::vector<std::string> termNamesColumn;
    std::vector<std::string> termIndicesColumn;

    int clusterNum = 1;
    for (const auto& cluster : finalClusters) {
        if (cluster.size() >= static_cast<size_t>(finalGroupMembership)) {
            std::vector<std::string> clusterTermNames;
            std::string termIndicesStr;

            for (int term_index : cluster) {
                clusterTermNames.push_back(Rcpp::as<std::string>(terms[term_index]));
                if (!termIndicesStr.empty()) {
                    termIndicesStr += ", ";
                }
                termIndicesStr += std::to_string(term_index);
            }

            clusterColumn.push_back(clusterNum++);
            termNamesColumn.push_back(StringUtils::vectorToString(clusterTermNames, ", "));
            termIndicesColumn.push_back(termIndicesStr);
        }
    }

    Rcpp::DataFrame clusters_df = Rcpp::DataFrame::create(
        Rcpp::Named("Cluster") = clusterColumn,
        Rcpp::Named("TermNames") = termNamesColumn,
        Rcpp::Named("TermIndices") = termIndicesColumn
    );

    return Rcpp::List::create(
        Rcpp::Named("clusters") = clusters_df
    );
}

void DavidClustering::calculateKappaScores() {
    for (int i = 0; i < n_terms; ++i) {
        for (int j = i + 1; j < n_terms; ++j) {
            std::unordered_set<std::string> term1_genes = StringUtils::splitStringToUnorderedSet(Rcpp::as<std::string>(geneIDs[i]), ",");
            std::unordered_set<std::string> term2_genes = StringUtils::splitStringToUnorderedSet(Rcpp::as<std::string>(geneIDs[j]), ",");

            int term1term2 = 0;
            for (const auto& gene : term1_genes) {
                if (term2_genes.count(gene)) {
                    term1term2++;
                }
            }

            int posTerm1Total = term1_genes.size();
            int posTerm2Total = term2_genes.size();
            int term1only = posTerm1Total - term1term2;
            int term2only = posTerm2Total - term1term2;
            int term1term2Non = totalGeneCount - term1term2 - term1only - term2only;

            double oab = static_cast<double>(term1term2 + term1term2Non) / totalGeneCount;
            double aab = (static_cast<double>(posTerm1Total) * posTerm2Total + static_cast<double>(totalGeneCount - posTerm1Total) * (totalGeneCount - posTerm2Total)) / (totalGeneCount * totalGeneCount);

            double kappa = (aab == 1) ? 1.0 : (oab - aab) / (1 - aab);

            kappaMatrix[i][j] = kappa;
            kappaMatrix[j][i] = kappa;
        }
    }
}

void DavidClustering::findInitialSeeds() {
    for (int i = 0; i < n_terms; ++i) {
        TermSet neighbors;
        for (int j = 0; j < n_terms; ++j) {
            if (i == j) continue;
            if (kappaMatrix[i][j] > similarityThreshold) {
                neighbors.insert(j);
            }
        }

        if (neighbors.size() >= static_cast<size_t>(initialGroupMembership - 1)) {
            TermSet current_seed = neighbors;
            current_seed.insert(i);

            int totalPairs = 0;
            int passedPair = 0;
            std::vector<int> seed_vec(current_seed.begin(), current_seed.end());
            for (size_t k = 0; k < seed_vec.size(); ++k) {
                for (size_t l = k + 1; l < seed_vec.size(); ++l) {
                    totalPairs++;
                    if (kappaMatrix[seed_vec[k]][seed_vec[l]] > similarityThreshold) {
                        passedPair++;
                    }
                }
            }

            if (totalPairs > 0 && (static_cast<double>(passedPair) / totalPairs) > multipleLinkageThreshold) {
                initialSeeds.push_back(current_seed);
            }
        }
    }
}

void DavidClustering::mergeSeeds() {
    std::list<TermSet> remainingSeeds(initialSeeds.begin(), initialSeeds.end());

    while (!remainingSeeds.empty()) {
        TermSet currentSeed = remainingSeeds.front();
        remainingSeeds.pop_front();

        while (true) {
            double bestScore = 0.0;
            auto bestIt = remainingSeeds.end();

            for (auto it = remainingSeeds.begin(); it != remainingSeeds.end(); ++it) {
                double score = calculateDiceCoefficient(currentSeed, *it);
                if (score > multipleLinkageThreshold && score > bestScore) {
                    bestScore = score;
                    bestIt = it;
                }
            }

            if (bestIt != remainingSeeds.end()) {
                currentSeed.insert(bestIt->begin(), bestIt->end());
                remainingSeeds.erase(bestIt);
            } else {
                break;
            }
        }
        finalClusters.push_back(currentSeed);
    }
}

double DavidClustering::calculateDiceCoefficient(const TermSet& seed1, const TermSet& seed2) {
    int commonCount = 0;
    for (int term : seed1) {
        if (seed2.count(term)) {
            commonCount++;
        }
    }
    return 2.0 * commonCount / (seed1.size() + seed2.size());
}

// Exported function to be called from R
// [[Rcpp::export]]
Rcpp::List runDavidClustering(
    Rcpp::CharacterVector terms,
    Rcpp::CharacterVector geneIDs,
    double similarityThreshold,
    int initialGroupMembership,
    int finalGroupMembership,
    double multipleLinkageThreshold) {

    DavidClustering david(
        terms,
        geneIDs,
        similarityThreshold,
        initialGroupMembership,
        finalGroupMembership,
        multipleLinkageThreshold
    );

    return david.run();
}
