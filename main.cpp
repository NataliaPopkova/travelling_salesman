#include "map.h"

int main(int argc, char** argv) {
    using namespace tsp;

    std::srand(std::time(nullptr));

    int N_cities = 10;
    Map map(N_cities);
    map.InitializeCities();

    std::ofstream file("Cities.txt");
    if (!file.is_open())
        std::cout << "Cities: File wasn't open!\n";
    for (int i = 0; i < map.GetWidth(); i++) {
        file << map.GetCities()[i].first << "    " << map.GetCities()[i].second
             << std::endl;
    }

    map.GetSeq() = Map::SimulatedAnnealing(map, 1000000000, 100);

    std::ofstream file1("Seq.txt");
    if (!file1.is_open())
        std::cout << "Seq: File wasn't open!\n";
    for (int i = 0; i < map.GetWidth(); i++) {
        file1 << map.GetSeq()[i] << std::endl;
    }
    file1 << map.GetSeq()[0] << std::endl;

    std::cout << "done!\n";
    return 0;
}
