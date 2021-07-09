/*
/   honeycombVentingSim.cpp
/   Calvin Boyle - 2021
/   Carnegie Mellon University
*/

#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
using namespace std;

# define M_PI           3.14159265358979323846  /* pi */

struct Cell
{
    double pressure;
    float temperature;
    double mass;
};

// Numerical Method Properties - currently euler's method - will update to RK4
float DT = 0.1;        // sec
float END_TIME = 100;     // sec
vector<vector<double>> data;
double max_pressure = 0;

// Falcon 9 Properties
double FAIRING_DEPRESSURE_RATE = 2800; // Pa per sec
double FAIRING_DEPRESSURE_STEP = FAIRING_DEPRESSURE_RATE * DT;

// Gas Properties
double R = 8.314462;        // (J/K-mol) universal gas constant
double R_AIR = 0.287058;    // (J/g-K) specific gas constant dry air
float M_AIR = .02897;       // molecular mass of dry air
float K_AIR = 1.4;          // isentropic exponent, ratio of specific heats Cp/Cv
double CRITICAL_PRESSURE_RATIO = (2 / pow((K_AIR + 1), (K_AIR / (K_AIR - 1))));    // choked flow transition point
float STARTING_TEMP = 300;  // K

// precomputes for speeeeeeeeed
float pre_1 = ((2 * K_AIR) / (K_AIR - 1)) * (M_AIR / R);
float pre_2 = 2 / K_AIR;
float pre_3 = (K_AIR + 1) / K_AIR;
float pre_4 = ((K_AIR * M_AIR) / R);
float pre_5 = pow((2 / (K_AIR + 1)), (K_AIR + 1) / (K_AIR - 1));

// Panel Properties
int NUM_ORIFACE = 8;
float ORIFACE_DIA = 0.002;         // m
float ORIFACE_DISCHARGE_COEFF = 0.63;   // current best guess TODO update to dynamic recalculation
double ORIFACE_AREA = (M_PI / 4.) * ORIFACE_DIA * ORIFACE_DIA;  // meter^2
double CHASSIS_VOLUME = 0.022;  // meter^3

double IdealGasPressure(float T, double m)
{
    // within 1% error for Air above -130 C
    // within 0.1% error at 23 C
    return (m * R_AIR * T) / CHASSIS_VOLUME;
}

double IdealGasMass(float T, double P)
{
    // within 1% error for Air above -130 C
    // within 0.1% error at 23 C
    return (P * CHASSIS_VOLUME) / (R_AIR * T);
}

double CalcMassFlow(double P_in, double P_out, double T)
{
    int flow_dir = 1;   // -1 into cell, 1 out of cell
    double temp;

    if (P_out > P_in)
    {
        flow_dir = -1;
        temp = P_out;
        P_out = P_in;
        P_in = temp;
    }

    double P_ratio = P_out / P_in;

    if (P_ratio > CRITICAL_PRESSURE_RATIO) // subsonic flow
        return (double)flow_dir * ORIFACE_DISCHARGE_COEFF * ORIFACE_AREA * P_in * pow(((pre_1 / T) * (pow(P_ratio, pre_2) - pow(P_ratio, pre_3))), 0.5);
    else // sonic (choked) flow
        return (double)flow_dir * ORIFACE_DISCHARGE_COEFF * ORIFACE_AREA * P_in * pow((pre_4 / T) * pre_5, 0.5);
}

void InitializeCell(Cell &Cell)
{
    // STP
    Cell.mass = IdealGasMass(300, 100000);
    Cell.pressure = 100000;
    Cell.temperature = STARTING_TEMP;
}

void EnforceBoundaryConditions(Cell &Cell)
{
    if (Cell.pressure > FAIRING_DEPRESSURE_STEP)
    {
        Cell.pressure -= FAIRING_DEPRESSURE_STEP;
    }
    else {
        Cell.pressure = 0;
    }
}

void SolveMassFlowPressure(Cell &Chassis, Cell &Fairing)
{
            Chassis.mass -= CalcMassFlow(Chassis.pressure, Fairing.pressure, Chassis.temperature) * NUM_ORIFACE;
            Chassis.pressure = IdealGasPressure(Chassis.temperature, Chassis.mass);
}

void WriteData(ofstream &f, float time, Cell &Chassis, Cell &Fairing)
{
    f << time << ", " << Fairing.pressure << ", " << Chassis.pressure << ", " << endl;
}

int main(void)
{
    /*
    1. initialize all cells
    Main Loop
    2. change all exterior cells by fairing rate
    3. loop through all interior cells
        3a. calculate mass flow in/out all 6 surrounding cells
        3b. sum mass flow rates
        3c. store mass flow rate in temp variable
    4. loop through all interior cells again
        4a. update cell mass by flow rate * dt
        4b. recalculate cell pressure
    */
    int VIEW_H = 1000;
    Cell Chassis;
    Cell Fairing;
    InitializeCell(Chassis);
    InitializeCell(Fairing);
    float time = 0;
    ofstream f("output.txt");
    int count = 0;

    while (time < END_TIME)
    {
        if (count == 0 || count == 50)
        {
            WriteData(f, time, Chassis, Fairing);
            cout << time << endl;
            count = 0;
        }
        count++;
        EnforceBoundaryConditions(Fairing);
        SolveMassFlowPressure(Chassis, Fairing);
        time += DT;
    }
    return 0;
}
