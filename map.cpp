#include "map.h"

Map::Map(size_t size) : size(size),
                        Cities(size, std::make_pair(0, 0)),
                        Seq(size)
{
}

size_t Map::GetWidth() const { return size; }  //возвращает количество столбцов
size_t Map::GetHeight() const { return size; } // возвращает количество строк

std::vector<std::pair<double, double>> &Map::GetCities() { return Cities; }
std::vector<int> &Map::GetSeq() { return Seq; }

void SwapSubsequence(int i, int j, std::vector<int> &arr)
{
    double a;
    if ((j - i - 1) % 2 != 0)
    { //если нечетное количество эл-тов между ними
        int z = i + (j - i - 2) / 2 + 1;
        for (int k = i; k <= j / 2; ++k)
        {
            a = arr[k];
            arr[k] = arr[k + 2 * (z - k)];
            arr[k + 2 * (z - k)] = a;
        }
    }
    else
    {
        int z = i + (j - i - 1) / 2;
        for (int k = i; k <= j / 2; ++k)
        {
            a = arr[k];
            arr[k] = arr[k + 2 * (z - k) + 1];
            arr[k + 2 * (z - k) + 1] = a;
        }
    }
}

void Map::InitializeCities()
{
    for (int i = 0; i < this->GetHeight(); ++i)
    {
        //{GetCities()[i] = {rand()%10, rand()%10};}
        {
            GetCities()[i] = {random(0, 100.), random(0, 100.)};
        }
    }
}

double Map::CalculateEnergy(std::vector<int> stateCandidate)
{
    int n = stateCandidate.size();
    double E = 0;
    for (int i = 1; i < n - 1; ++i)
    {
        E = E + sqrt((GetCities()[i + 1].first - GetCities()[i].first) * (GetCities()[i + 1].first - GetCities()[i].first) + (GetCities()[i + 1].second - GetCities()[i].second) * (GetCities()[i + 1].second - GetCities()[i].second));
    }
    E = E + sqrt((GetCities()[n - 1].first - GetCities()[n - 2].first) * (GetCities()[n - 1].first - GetCities()[n - 2].first) + (GetCities()[n - 1].second - GetCities()[n - 2].second) * (GetCities()[n - 1].second - GetCities()[n - 2].second));
    return E;
}

void Map::GenerateStateCandidate(std::vector<int> &)
{
    //    seq - предыдущее состояние (маршрут), из которого мы хотим получить состояние-кандидат
    int n = GetWidth(); // определяем размер последовательности
    int i = rand() % n; //генерируем целое случайное число
    int j = rand() % n; // генерируем целое случайное число

    if (i > j)
        SwapSubsequence(j, i, GetSeq()); // обращаем подпоследовательность
    else
        SwapSubsequence(i, j, GetSeq()); // обращаем подпоследовательность
    //return GetSeq();
}

double DecreaseTemperature(double initialTemperature, int i)
{
    //initialTemperature - начальная температура
    //i - номер итерации
    return initialTemperature * 0.1 / i;
}

double GetTransitionProbability(double dE, double T)
{
    return exp(-dE / T);
}

double MakeTransit(double probability)
{
    double value = std::rand() % 1;
    if (value <= probability)
        return 1;
    else
        return 0;
}

double random(double min, double max)
{
    return (double)(rand()) / RAND_MAX * (max - min) + min;
}

std::vector<int> SimulatedAnnealing(Map map, double initialTemperature, double endTemperature)
{

    int n = map.GetHeight(); // получаем размер вектора городов
                             //    for (int i = 0; i < n; ++i) {
                             //        map.GetSeq()[i] = (rand() % 100);
                             //    }

    for (int i = 0; i < n; ++i)
    {
        map.GetSeq()[i] = i;
    }

    for (int i = 0; i < n; ++i)
    {
        //map.GetSeq()[i] = (rand() % 100);
        random_shuffle(map.GetSeq().begin(), map.GetSeq().end());
    }
    // Функция randperm(n) - генерирует случайныую последовательность из целых чисел от 1 до n

    double currentEnergy = map.CalculateEnergy(map.GetSeq()); // вычисляем энергию для первого состояния
    double T = initialTemperature;

    for (int i = 1; i < 1000000; ++i)
    {
        map.GenerateStateCandidate(map.GetSeq()); // получаем состояние-кандидат
        std::vector<int> stateCandidate = map.GetSeq();
        double candidateEnergy = map.CalculateEnergy(stateCandidate); // вычисляем его энергию

        double p = GetTransitionProbability(candidateEnergy - currentEnergy, T);
        if (candidateEnergy < currentEnergy)
        {                                    // если кандидат обладает меньшей энергией
            currentEnergy = candidateEnergy; // то оно становится текущим состоянием
            map.GetSeq() = stateCandidate;
        }
        else
        { // иначе, считаем вероятность
            if (MakeTransit(p))
            { // и смотрим, осуществится ли переход
                currentEnergy = candidateEnergy;
                map.GetSeq() = stateCandidate;
            }
        }

        T = DecreaseTemperature(initialTemperature, i); // уменьшаем температуру
        if (T <= endTemperature)                        // условие выхода
            return map.GetSeq();
    }

    //return map.GetSeq();
}
