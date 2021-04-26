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
float DT = 0.01;        // sec
float END_TIME = 10;     // sec
vector<vector<double>> data;
double max_pressure = 0;

// Falcon 9 Properties
float FAIRING_DEPRESSURE_RATE = 2800; // Pa per sec
float FAIRING_DEPRESSURE_STEP = FAIRING_DEPRESSURE_RATE * DT;

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
float ORIFACE_DIA = 0.0005;         // m
float ORIFACE_DISCHARGE_COEFF = 0.63;   // current best guess TODO update to dynamic recalculation
double ORIFACE_AREA = (M_PI / 4.) * ORIFACE_DIA * ORIFACE_DIA;  // meter^2
float PANEL_WIDTH = 600;                // mm
float PANEL_HEIGHT = 600;               // mm
float CORE_THICK = 10;                  // mm
float CELL_SIZE = 6.35;                 // mm
double CELL_HEX_AREA = 6. * (CELL_SIZE / 2000.) * (CELL_SIZE / 2000.) * tan(M_PI / 6.); // meter^2
double CELL_VOLUME = (CORE_THICK / 1000.) * CELL_HEX_AREA;  // meter^3
int num_cell_x = PANEL_WIDTH / CELL_SIZE + 2;
int num_cell_y = PANEL_HEIGHT / CELL_SIZE + 2;

double IdealGasPressure(float T, double m)
{
    // within 1% error for Air above -130 C
    // within 0.1% error at 23 C
    return (m * R_AIR * T) / CELL_VOLUME;
}

double IdealGasMass(float T, double P, float V)
{
    // within 1% error for Air above -130 C
    // within 0.1% error at 23 C
    return (P * V) / (R_AIR * T);
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

void InitializePanel(vector<vector<Cell>> &Panel)
{
    vector<Cell> newRow;

    for (int i = 0; i < num_cell_y; i++)
    {
        for (int j = 0; j < num_cell_x; j++)
        {
            Cell newCell;
            newCell.mass = IdealGasMass(300, 100000, CELL_VOLUME);
            newCell.pressure = 100000;
            newCell.temperature = STARTING_TEMP;
            newRow.push_back(newCell);
        }
        Panel.push_back(newRow);
        newRow.clear();
    }
}

void InitializePanel(vector<vector<double>>& Panel)
{
    vector<double> newRow;

    for (int i = 0; i < num_cell_y; i++)
    {
        for (int j = 0; j < num_cell_x; j++)
        {
            double mass = 0;
            newRow.push_back(mass);
        }
        Panel.push_back(newRow);
        newRow.clear();
    }
}

void EnforceBoundaryConditions(vector<vector<Cell>> &Panel)
{
    if (Panel[0][0].pressure > FAIRING_DEPRESSURE_STEP)
    {
        // left-right edges
        for (int i = 0; i < num_cell_y; i++)
        {
            Panel[i][0].pressure -= (double)FAIRING_DEPRESSURE_STEP;
            Panel[i][num_cell_x - 1].pressure -= (double)FAIRING_DEPRESSURE_STEP;
        }

        // top-bottom edges
        for (int j = 1; j < num_cell_x-1; j++)
        {
            Panel[0][j].pressure -= (double)FAIRING_DEPRESSURE_STEP;
            Panel[num_cell_y - 1][j].pressure -= (double)FAIRING_DEPRESSURE_STEP;
        }
    }
}

void SolveMassFlow(vector<vector<Cell>>& Panel, vector<vector<double>>& Panel_delta)
{
    // only for interior cells
    for (int i = 1; i < num_cell_y - 1; i++)
    {
        for (int j = 1; j < num_cell_x - 1; j++)
        {
            // sum mass flow for all 6 neighbor cells - not optimized b/c double calcs every oriface in both directions
            Panel_delta[i][j] = 0.0;
            Panel_delta[i][j] += CalcMassFlow(Panel[i][j].pressure, Panel[i-1][j-1].pressure, Panel[i][j].temperature);
            Panel_delta[i][j] += CalcMassFlow(Panel[i][j].pressure, Panel[i-1][j].pressure, Panel[i][j].temperature);
            Panel_delta[i][j] += CalcMassFlow(Panel[i][j].pressure, Panel[i][j+1].pressure, Panel[i][j].temperature);
            Panel_delta[i][j] += CalcMassFlow(Panel[i][j].pressure, Panel[i+1][j].pressure, Panel[i][j].temperature);
            Panel_delta[i][j] += CalcMassFlow(Panel[i][j].pressure, Panel[i+1][j-1].pressure, Panel[i][j].temperature);
            Panel_delta[i][j] += CalcMassFlow(Panel[i][j].pressure, Panel[i][j-1].pressure, Panel[i][j].temperature);
        }
    }
}

void UpdateMassPressure(vector<vector<Cell>>& Panel, vector<vector<double>>& Panel_delta)
{
    max_pressure = 0;
    for (int i = 1; i < num_cell_y - 1; i++)
    {
        for (int j = 1; j < num_cell_x - 1; j++)
        {
            Panel[i][j].mass -= Panel_delta[i][j] * DT;
            Panel[i][j].pressure = IdealGasPressure(Panel[i][j].temperature, Panel[i][j].mass);

            if (Panel[i][j].pressure > max_pressure)
                max_pressure = Panel[i][j].pressure;
        }
    }
}

void WriteData(ofstream &f, float time, vector<vector<Cell>> &Panel)
{
    f << time << ", " << Panel[0][0].pressure << ", " << Panel[1][1].pressure << ", " << max_pressure << endl;
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
    vector<vector<Cell>> Panel;
    vector<vector<double>> Panel_delta;
    InitializePanel(Panel);
    InitializePanel(Panel_delta);
    float time = 0;
    ofstream f("output.txt");
    int count = 0;

    FsOpenWindow(16, 16, VIEW_H, VIEW_H, 0, "Test");
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glColor3ub(100, 100, 100);
    glPointSize(30);

    while (time < END_TIME)
    {
        if (count == 0 || count == 50)
        {
            WriteData(f, time, Panel);
            cout << time << endl;
            count = 0;
        }
        count++;
        EnforceBoundaryConditions(Panel);
        SolveMassFlow(Panel, Panel_delta);
        UpdateMassPressure(Panel, Panel_delta);
        time += DT;
    }
    return 0;
}
