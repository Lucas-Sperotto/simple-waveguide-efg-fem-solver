#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>
#include <algorithm>

const double pi = 3.14159265358979323846;
const double c = 299792458.0; // m/s

struct Modo {
    std::string tipo; // "TM" ou "TE"
    int m;
    int n;
    double raiz;
    double fc_GHz;
};

int main() {
    std::ifstream infile("bessel_roots.csv");
    if (!infile) {
        std::cerr << "Erro ao abrir o arquivo bessel_roots.csv.\n";
        return 1;
    }

    double a = 0.01; // Raio do guia de onda (m)
    std::vector<Modo> modos;
    std::string line;

    // Ignora a primeira linha (cabeçalho)
    std::getline(infile, line);

    while (std::getline(infile, line)) {
        std::istringstream ss(line);
        std::string tipo_str, m_str, n_str, raiz_str;

        std::getline(ss, tipo_str, ',');
        std::getline(ss, m_str, ',');
        std::getline(ss, n_str, ',');
        std::getline(ss, raiz_str, ',');

        Modo modo;
        modo.tipo = tipo_str;
        modo.m = std::stoi(m_str);
        modo.n = std::stoi(n_str);
        modo.raiz = std::stod(raiz_str);

        double fc_Hz = (modo.raiz * c) / (2 * pi * a);
        modo.fc_GHz = fc_Hz / 1e9;

        modos.push_back(modo);
    }

    infile.close();

    // Ordena por frequência de corte
    std::sort(modos.begin(), modos.end(), [](const Modo& a, const Modo& b) {
        return a.fc_GHz < b.fc_GHz;
    });

    std::ofstream outfile("frequencias_de_corte.csv");
    if (!outfile) {
        std::cerr << "Erro ao criar o arquivo frequencias_de_corte.csv.\n";
        return 1;
    }

    outfile << std::fixed << std::setprecision(15);
    outfile << "Modo,m,n,Frequencia_de_corte_GHz\n";

    for (const auto& modo : modos) {
        outfile << modo.tipo << "," << modo.m << "," << modo.n << "," << modo.fc_GHz << "\n";
    }

    outfile.close();

    std::cout << "Arquivo 'frequencias_de_corte.csv' gerado com sucesso!\n";
    return 0;
}
