#include "map.h"

namespace tsp {

Map::Map(size_t size)
        : size(size), Cities(size, std::make_pair(0, 0)), Seq(size) {}

size_t Map::GetWidth() const {
    return size;
}  //возвращает количество столбцов
size_t Map::GetHeight() const {
    return size;
}  // возвращает количество строк

std::vector<std::pair<double, double>>& Map::GetCities() {
    return Cities;
}
std::vector<int>& Map::GetSeq() {
    return Seq;
}

void Map::InitializeCities() {
    for (int i = 0; i < this->GetHeight(); ++i) {
        { GetCities()[i] = {random(0, 100.), random(0, 100.)}; }
    }
}

double Map::CalculateEnergy(std::vector<int> stateCandidate) {
    int    n = stateCandidate.size();
    double E = 0;
    for (int i = 1; i < n - 1; i++) {
        E += sqrt(pow((GetCities()[stateCandidate[i + 1]].first -
                       GetCities()[stateCandidate[i]].first),
                      2) +
                  pow((GetCities()[stateCandidate[i + 1]].second -
                       GetCities()[stateCandidate[i]].second),
                      2));
    }

    E += sqrt(pow((GetCities()[stateCandidate[n - 1]].first -
                   GetCities()[stateCandidate[n - 2]].first),
                  2) +
              pow((GetCities()[stateCandidate[n - 1]].second -
                   GetCities()[stateCandidate[n - 2]].second),
                  2));

    return E;
}

void Map::GenerateStateCandidate(std::vector<int>& seq) {
    //    seq - предыдущее состояние (маршрут), из которого мы хотим получить
    //    состояние-кандидат
    int n = seq.size();  // определяем размер последовательности
    int i = rand() % n;  //генерируем целое случайное число
    int j = rand() % n;  // генерируем целое случайное число

    auto seq_iterator_i = seq.begin();
    auto seq_iterator_j = seq.begin();

    std::advance(seq_iterator_i, i);
    std::advance(seq_iterator_j, j);

    if (i > j)
        std::reverse(seq_iterator_j,
                     seq_iterator_i);  // обращаем подпоследовательность
    else
        std::reverse(seq_iterator_i,
                     seq_iterator_j);  // обращаем подпоследовательность
}

double Map::DecreaseTemperature(double& initialTemperature, int i) {
    // initialTemperature - начальная температура
    // i - номер итерации
    return initialTemperature * 0.1 / i;
}

double Map::GetTransitionProbability(double dE, double T) {
    return exp(-dE / T);
}

bool Map::MakeTransit(double probability) {
    double value = std::rand() % 1;
    if (value <= probability)
        return true;
    else
        return false;
}

double Map::random(double min, double max) {
    return (double)(rand()) / RAND_MAX * (max - min) + min;
}

std::vector<int> Map::SimulatedAnnealing(Map map, double initialTemperature,
                                         double endTemperature) {
    int n = map.GetHeight();  // получаем размер вектора городов

    for (int i = 0; i < n; ++i) {
        map.GetSeq()[i] = i;
    }

    // задаём начальное состояние, как случайный маршрут
    std::random_device rd;
    std::mt19937       g(rd());
    for (int i = 0; i < n; ++i) {
        std::shuffle(map.GetSeq().begin(), map.GetSeq().end(), g);
    }

    double currentEnergy = map.CalculateEnergy(
        map.GetSeq());  // вычисляем энергию для первого состояния

    double T = initialTemperature;

    // for (int i = 1; i < 10000; ++i)
    int i = 0;
    while (T != endTemperature) {
        map.GenerateStateCandidate(
            map.GetSeq());  // получаем состояние-кандидат
        std::vector<int> stateCandidate(map.GetSeq());
        double           candidateEnergy =
            map.CalculateEnergy(stateCandidate);  // вычисляем его энергию

        double p = GetTransitionProbability(candidateEnergy - currentEnergy, T);
        if (candidateEnergy <
            currentEnergy) {  // если кандидат обладает меньшей энергией
            currentEnergy =
                candidateEnergy;  // то оно становится текущим состоянием
            map.GetSeq().assign(stateCandidate.begin(), stateCandidate.end());
        } else {  // иначе, считаем вероятность
            if (MakeTransit(p)) {  // и смотрим, осуществится ли переход
                currentEnergy = candidateEnergy;
                map.GetSeq().assign(stateCandidate.begin(),
                                    stateCandidate.end());
            }
        }

        i++;
        T = DecreaseTemperature(initialTemperature,
                                i);  // уменьшаем температуру
        if (T <= endTemperature)     // условие выхода
            return map.GetSeq();
    }

    return map.GetSeq();
}


}  // namespace tsp