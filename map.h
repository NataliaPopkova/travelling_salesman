#include <math.h>
#include <algorithm>
#include <ctime>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

namespace tsp {

class Map {
public:
    Map() = default;
    Map(size_t size);

    size_t GetWidth() const;  //возвращает количество столбцов
    size_t GetHeight() const;  // возвращает количество строк

    std::vector<std::pair<double, double>>& GetCities();
    std::vector<int>&                       GetSeq();

    void   InitializeCities();
    double CalculateEnergy(std::vector<int> stateCandidate);
    void   GenerateStateCandidate(std::vector<int>&);

    static double DecreaseTemperature(double& initialTemperature, int i);
    static double GetTransitionProbability(double dE, double T);
    static double MakeTransit(double probability);
    static double random(double min, double max);
    static std::vector<int> SimulatedAnnealing(Map    map,
                                               double initialTemperature,
                                               double endTemperature);

private:
    std::vector<std::pair<double, double>> Cities;
    std::vector<int>                       Seq;
    size_t                                 size;
};

}  // namespace tsp
