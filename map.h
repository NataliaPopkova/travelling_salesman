#ifndef TSP_H
#define TSP_H

#include <random>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <algorithm>
#include <ctime>

void SwapSubsequence(int i, int j, std::vector<double> &arr);
double DecreaseTemperature(double initialTemperature, int i);
double GetTransitionProbability(double dE, double T);
double MakeTransit(double probability);
double random(double min, double max);

class Map
{

public:
    Map() = default;
    Map(size_t size);
    size_t GetWidth() const;  //возвращает количество столбцов
    size_t GetHeight() const; // возвращает количество строк

    std::vector<std::pair<double, double>> &GetCities();
    std::vector<int> &GetSeq();

    void InitializeCities();
    double CalculateEnergy(std::vector<int> stateCandidate);
    void GenerateStateCandidate(std::vector<int> &);

private:
    std::vector<std::pair<double, double>> Cities;
    std::vector<int> Seq;
    size_t size;
};
std::vector<int> SimulatedAnnealing(Map map, double initialTemperature, double endTemperature);
#endif
