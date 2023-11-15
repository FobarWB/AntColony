#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <random>

struct Ant {
    std::vector<int> Trail;
    std::vector<bool> visited;
    std::mt19937 random;

    Ant(int numberOfCities) : Trail(numberOfCities), visited(numberOfCities, false), random(std::random_device{}()) {}

    void Clear() {
        std::fill(Trail.begin(), Trail.end(), -1);
        std::fill(visited.begin(), visited.end(), false);
    }

    void VisitCity(int currentIndex, int city) {
        if (currentIndex + 1 < Trail.size() && city >= 0 && city < visited.size()) {
            Trail[currentIndex + 1] = city;
            visited[city] = true;
        }
    }

    bool Visited(int i) const {
        return visited[i];
    }

    double TrailLength(const std::vector<std::vector<double>>& graph) const {
        double length = 0.0;
        for (size_t i = 0; i < Trail.size() - 1; i++) {
            length += graph[Trail[i]][Trail[i + 1]];
        }
        length += graph[Trail.back()][Trail[0]];
        return length;
    }
};

struct AntColony {
    std::vector<int> bestTour;
    std::vector<std::vector<double>> graph;
    std::vector<std::vector<double>> trails;
    std::vector<double> probabilities;
    std::vector<Ant> ants;
    int numberOfCities;
    int numberOfAnts;
    double evaporation;
    double alpha;
    double beta;
    double Q;
    std::mt19937 random;
    int currentIndex;
    std::vector<int> bestTourOrder;
    double bestTourLength;

    AntColony(std::vector<std::vector<double>>& graph, int numberOfAnts, double evaporation, double alpha, double beta, double Q, std::vector<int>& bestTour)
        : graph(graph), numberOfCities(graph.size()), numberOfAnts(numberOfAnts), evaporation(evaporation), alpha(alpha), beta(beta), Q(Q),
        random(std::random_device{}()), currentIndex(0), bestTourOrder(), bestTourLength(std::numeric_limits<double>::max()) {

        trails.resize(numberOfCities, std::vector<double>(numberOfCities, 1.0));
        probabilities.resize(numberOfCities);
        ants.resize(numberOfAnts, Ant(numberOfCities));
        this->bestTour = bestTour;
    }

    std::vector<int> CalculateAOC(int maxIterations) {
        std::vector<int> globalBestTour;
        double globalBestLength = std::numeric_limits<double>::max();

        for (int iteration = 0; iteration < maxIterations; iteration++) {
            SetupAnts();
            MoveAnts();
            UpdateTrails();
            UpdateBest();

            if (bestTourLength < globalBestLength) {
                globalBestLength = bestTourLength;
                globalBestTour = bestTourOrder;
            }
            std::cout<<"-------------------------------" << std::endl;
            std::cout << "Iteration " << iteration + 1 << ": Best Tour Length: " << bestTourLength << std::endl;
            std::cout << "Best Tour: ";
            for (int cityIndex : bestTourOrder) {
                std::cout << cityIndex + 1 << " - ";
            }
            std::cout << bestTourOrder[0]+1;
            std::cout << std::endl << std::endl;
        }
        return globalBestTour;
    }

private:
    void SetupAnts() {
        for (auto& ant : ants) {
            ant.Clear();
            ant.VisitCity(-1, random() % numberOfCities);
        }
        currentIndex = 0;
    }

    void MoveAnts() {
        for (int i = currentIndex; i < numberOfCities - 1; i++) {
            for (auto& ant : ants) {
                ant.VisitCity(currentIndex, SelectNextCity(ant));
            }
            currentIndex++;
        }
    }

    int SelectNextCity(Ant& ant) {
        int i = ant.Trail[currentIndex];
        double pheromone = 0.0;

        for (int l = 0; l < numberOfCities; l++) {
            if (!ant.Visited(l)) {
                if (i >= 0 && i < numberOfCities && l >= 0 && l < numberOfCities) {
                    pheromone += std::pow(trails[i][l], alpha) * std::pow(1.0 / graph[i][l], beta);
                }
            }
        }

        for (int j = 0; j < numberOfCities; j++) {
            if (ant.Visited(j)) {
                probabilities[j] = 0.0;
            }
            else {
                if (i >= 0 && i < numberOfCities && j >= 0 && j < numberOfCities) {
                    double numerator = std::pow(trails[i][j], alpha) * std::pow(1.0 / graph[i][j], beta);
                    probabilities[j] = numerator / pheromone;
                }
            }
        }

        double r = static_cast<double>(random()) / random.max();
        double total = 0;

        for (int q = 0; q < numberOfCities; q++) {
            total += probabilities[q];
            if (total >= r) {
                return q;
            }
        }

        // This code should never be reached, but to avoid compiler errors
        return -1;
    }

    void UpdateTrails() {
        for (int i = 0; i < numberOfCities; i++) {
            for (int j = 0; j < numberOfCities; j++) {
                if (i >= 0 && i < numberOfCities && j >= 0 && j < numberOfCities) {
                    trails[i][j] *= evaporation;
                }
                else {
                    std::cerr << "Error: Indices out of bounds in trails[" << i << "][" << j << "]" << std::endl;
                }
            }
        }

        for (const auto& ant : ants) {
            double contribution = Q / ant.TrailLength(graph);
            for (size_t i = 0; i < ant.Trail.size() - 1; i++) {
                if (ant.Trail[i] >= 0 && ant.Trail[i] < numberOfCities && ant.Trail[i + 1] >= 0 && ant.Trail[i + 1] < numberOfCities) {
                    trails[ant.Trail[i]][ant.Trail[i + 1]] += contribution;
                }
                else {
                    std::cerr << "Error: Indices out of bounds in trails[" << ant.Trail[i] << "][" << ant.Trail[i + 1] << "]" << std::endl;
                }
            }
            if (ant.Trail.back() >= 0 && ant.Trail.back() < numberOfCities && ant.Trail[0] >= 0 && ant.Trail[0] < numberOfCities) {
                trails[ant.Trail.back()][ant.Trail[0]] += contribution;
            }
            else {
                std::cerr << "Error: Indices out of bounds in trails[" << ant.Trail.back() << "][" << ant.Trail[0] << "]" << std::endl;
            }
        }
    }

    void UpdateBest() {
        if (bestTourOrder.empty()) {
            bestTourOrder = ants[0].Trail;
            bestTourLength = ants[0].TrailLength(graph);
        }

        for (const auto& ant : ants) {
            if (ant.TrailLength(graph) < bestTourLength) {
                bestTourLength = ant.TrailLength(graph);
                bestTourOrder = ant.Trail;
            }
        }
    }
};

double CalculateTourLength(const std::vector<int>& tour, const std::vector<std::vector<double>>& distances) {
    double length = 0.0;
    for (size_t i = 0; i < tour.size() - 1; i++) {
        length += distances[tour[i]][tour[i + 1]];
    }
    length += distances[tour.back()][tour[0]];  // Add the distance from the last city to the starting one
    return length;
}

int main() {
    std::vector<std::vector<double>> distances = {
        {0, 3, 5, 6},
        {6, 0, 9, 3},
        {5, 9, 0, 12},
        {6, 3, 12, 0}
    };

    int numAnts = 5;
    double evaporation = 0.1;
    double alpha = 1.0;
    double beta = 1.0;
    double Q = 20;
    std::vector<int> bestTour;
    AntColony antColony(distances, numAnts, evaporation, alpha, beta, Q, bestTour);

    std::vector<int> bestTourResult = antColony.CalculateAOC(numAnts);
    double bestTourLength = CalculateTourLength(bestTourResult, distances);
    return 0;
}