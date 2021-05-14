#include "map.h"

namespace tsp
{

    Map::Map(size_t size)
        : size(size), Cities(size, std::make_pair(0, 0)), Seq(size) {}

    size_t Map::GetWidth() const
    {
        return size;
    } //возвращает количество столбцов
    size_t Map::GetHeight() const
    {
        return size;
    } // возвращает количество строк

    std::vector<std::pair<double, double>> &Map::GetCities()
    {
        return Cities;
    }
    std::vector<int> &Map::GetSeq()
    {
        return Seq;
    }

    void Map::SetSeq(std::vector<int> candidate)
    {
        Seq = candidate;
    };

    void Map::InitializeCities()
    {
        for (int i = 0; i < this->GetHeight(); ++i)
        {
            {
                GetCities()[i] = {random(0, 100.), random(0, 100.)};
            }
        }
    }

    double Map::CalculateEnergy(std::vector<int> stateCandidate)
    {
        int n = stateCandidate.size();
        double E = 0;
        for (int i = 0; i < n - 1; i++)
        {
            E += sqrt(pow((GetCities()[stateCandidate[i + 1]].first -
                           GetCities()[stateCandidate[i]].first),
                          2) +
                      pow((GetCities()[stateCandidate[i + 1]].second -
                           GetCities()[stateCandidate[i]].second),
                          2));
        }

        E += sqrt(pow((GetCities()[stateCandidate[n - 1]].first -
                       GetCities()[stateCandidate[0]].first),
                      2) +
                  pow((GetCities()[stateCandidate[n - 1]].second -
                       GetCities()[stateCandidate[0]].second),
                      2));

        return E;
    }

    std::vector<int> Map::GenerateStateCandidate(std::vector<int> &candidate)
    {
        //    seq - предыдущее состояние (маршрут), из которого мы хотим получить
        //    состояние-кандидат
        std::vector<int> seq{candidate};
        int n = seq.size(); // определяем размер последовательности
        int i = rand() % n; // генерируем целое случайное число
        int j = rand() % n; // генерируем целое случайное число

        auto seq_iterator_i = seq.begin();
        auto seq_iterator_j = seq.begin();

        std::advance(seq_iterator_i, i);
        std::advance(seq_iterator_j, j);

        if (i > j)
            std::reverse(seq_iterator_j,
                         seq_iterator_i); // обращаем подпоследовательность
        else
            std::reverse(seq_iterator_i,
                         seq_iterator_j); // обращаем подпоследовательность

        return seq;
    }

    double Map::DecreaseTemperature(double &initialTemperature, int i)
    {
        // initialTemperature - начальная температура
        // i - номер итерации
        return initialTemperature * 0.1 / i;
    }

    double Map::GetTransitionProbability(double dE, double T)
    {
        return expf(-dE / T);
    }

    bool Map::MakeTransit(double probability)
    {
        double value = (float)rand() / RAND_MAX;
        if (value <= probability)
            return true;
        else
            return false;
    }

    double Map::random(double min, double max)
    {
        return (double)(rand()) / RAND_MAX * (max - min) + min;
    }

    std::vector<int> Map::SimulatedAnnealing(Map map, double initialTemperature,
                                             double endTemperature)
    {
        int n = map.GetHeight(); // получаем размер вектора городов

        std::vector<int> seq(n);
        for (int i = 0; i < n; ++i)
        {
            seq[i] = i;
        }
        map.SetSeq(seq);

        // задаём начальное состояние, как случайный маршрут
        std::random_device rd;
        std::mt19937 g(rd());
        for (int i = 0; i < n; ++i)
        {
            std::shuffle(map.GetSeq().begin(), map.GetSeq().end(), g);
        }

        double currentEnergy = map.CalculateEnergy(map.GetSeq());
        // вычисляем энергию для первого состояния

        double T = initialTemperature;

        // for (int i = 1; i < 10000; ++i)

        int norm_transit = 0;
        int prob_transit = 0;
        int no_transit = 0;

        int i = 0;
        double p = 0.0;
        while (T != endTemperature)
        {
            // получаем состояние-кандидат
            std::vector<int> stateCandidate = map.GenerateStateCandidate(
                map.GetSeq());
            double candidateEnergy =
                map.CalculateEnergy(stateCandidate); // вычисляем его энергию

            if (candidateEnergy <
                currentEnergy)
            { // если кандидат обладает меньшей энергией
                currentEnergy =
                    candidateEnergy; // то оно становится текущим состоянием
                map.SetSeq(stateCandidate);
                norm_transit++;
            }
            else
            { // иначе, считаем вероятность
                p = GetTransitionProbability(candidateEnergy, T);

                if (MakeTransit(p))
                { // и смотрим, осуществится ли переход
                    currentEnergy = candidateEnergy;
                    map.SetSeq(stateCandidate);
                    prob_transit++;
                }
                else
                    no_transit++;
            }

            i++;
            T = DecreaseTemperature(initialTemperature,
                                    i); // уменьшаем температуру

            if (i == 1)
            {
                std::cout << "E = " << map.CalculateEnergy(map.GetSeq()) << "\n";
            }
        }

        std::cout << "E = " << map.CalculateEnergy(map.GetSeq()) << "\n";

        std::cout << "norm_transit = " << norm_transit << "\n";
        std::cout << "prob_transit = " << prob_transit << "\n";
        std::cout << "no_transit = " << no_transit << "\n";

        return map.GetSeq();
    }

} // namespace tsp