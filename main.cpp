#include "map.h"

int main(int argc, char** argv) {
    std::srand(std::time(nullptr));
    Map map(100);
    map.InitializeCities();

    std::ofstream file("Cities.txt");
    if (!file.is_open())
        std::cout << "Cities: File wasn't open!\n";
    for (int i = 0; i < map.GetWidth(); ++i) {
        file << map.GetCities()[i].first << "    " << map.GetCities()[i].second
             << std::endl;
    }

    map.GetSeq() = SimulatedAnnealing(map, 10, 0.00001);

    std::ofstream file1("Seq.txt");
    if (!file1.is_open())
        std::cout << "Seq: File wasn't open!\n";
    for (int i = 0; i < map.GetWidth(); ++i) {
        file1 << map.GetSeq()[i] << std::endl;
    }
    file1 << map.GetSeq()[0] << std::endl;

    // for (int i = 0; i < map.GetHeight(); ++i)
    // {
    //     std::cout << map.GetSeq()[i] << std::endl;
    // }

    std::cout << "done!";
    return 0;
}
