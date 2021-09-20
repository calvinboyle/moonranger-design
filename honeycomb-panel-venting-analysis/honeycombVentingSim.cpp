/*
/   honeycombVentingSim.cpp
/   Calvin Boyle - 2021
/   Carnegie Mellon University
*/

#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include "fssimplewindow.h"
#include "DrawingUtilNG.h"
using namespace std;

# define M_PI           3.14159265358979323846  /* pi */

struct Cell
{
    double pressure;
    float temperature;
    double mass;
};

// Numerical Method Properties - currently euler's method - will update to RK4
float dt = 0.01;        // sec
float end_time = 10;     // sec
vector<vector<double>> data;
double max_pressure = 0;

// Falcon 9 Properties
float fairing_depressure_rate = 2800; // Pa per sec
float fairing_depressure_increment = fairing_depressure_rate * dt;

// Gas Properties
double R = 8.314462;     // (J/K-mol) universal gas constant
double R_air = 0.287058; // (J/g-K) specific gas constant dry air
float M_air = .02897;   // molecular mass of dry air
float k_air = 1.4;      // isentropic exponent, ratio of specific heats Cp/Cv
double critical_pressure_ratio = (2 / pow((k_air + 1), (k_air / (k_air - 1))));    // choked flow transition point
float starting_temp = 300;  // K

// precomputes for speeeeeeeeed
float pre_1 = ((2 * k_air) / (k_air - 1)) * (M_air / R);
float pre_2 = 2 / k_air;
float pre_3 = (k_air + 1) / k_air;
float pre_4 = ((k_air * M_air) / R);
float pre_5 = pow((2 / (k_air + 1)), (k_air + 1) / (k_air - 1));

// Panel Properties
float oriface_diameter = 0.0005;         // m
float oriface_discharge_coeff = 0.63;   // current best guess TODO update to dynamic recalculation
double oriface_area = (M_PI / 4.) * oriface_diameter * oriface_diameter;  // meter^2
float panel_width = 600;                // mm
float panel_height = 600;               // mm
float core_thick = 10;                  // mm
float cell_size = 6.35;                 // mm
double cell_hex_area = 6. * (cell_size / 2000.) * (cell_size / 2000.) * tan(M_PI / 6.); // meter^2
double cell_volume = (core_thick / 1000.) * cell_hex_area;  // meter^3
int num_cell_x = panel_width / cell_size + 2;
int num_cell_y = panel_height / cell_size + 2;

double IdealGasPressure(float T, double m)
{
    // within 1% error for Air above -130 C
    // within 0.1% error at 23 C
    return (m * R_air * T) / cell_volume;
}

double IdealGasMass(float T, double P, float V)
{
    // within 1% error for Air above -130 C
    // within 0.1% error at 23 C
    return (P * V) / (R_air * T);
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

    if (P_ratio > critical_pressure_ratio) // subsonic flow
        return (double)flow_dir * oriface_discharge_coeff * oriface_area * P_in * pow(((pre_1 / T) * (pow(P_ratio, pre_2) - pow(P_ratio, pre_3))), 0.5);
    else // sonic (choked) flow
        return (double)flow_dir * oriface_discharge_coeff * oriface_area * P_in * pow((pre_4 / T) * pre_5, 0.5);
}

void InitializePanel(vector<vector<Cell>> &Panel)
{
    vector<Cell> newRow;

    for (int i = 0; i < num_cell_y; i++)
    {
        for (int j = 0; j < num_cell_x; j++)
        {
            Cell newCell;
            newCell.mass = IdealGasMass(300, 100000, cell_volume);
            newCell.pressure = 100000;
            newCell.temperature = starting_temp;
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
    if (Panel[0][0].pressure > fairing_depressure_increment)
    {
        // left-right edges
        for (int i = 0; i < num_cell_y; i++)
        {
            Panel[i][0].pressure -= (double)fairing_depressure_increment;
            Panel[i][num_cell_x - 1].pressure -= (double)fairing_depressure_increment;
        }

        // top-bottom edges
        for (int j = 1; j < num_cell_x-1; j++)
        {
            Panel[0][j].pressure -= (double)fairing_depressure_increment;
            Panel[num_cell_y - 1][j].pressure -= (double)fairing_depressure_increment;
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
            Panel[i][j].mass -= Panel_delta[i][j] * dt;
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

void Render(vector<vector<Cell>> Panel)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glBegin(GL_POINTS);

    double red, green, blue;

    for (int i = 0; i < num_cell_y; i++)
    {
        for (int j = 0; j < num_cell_x; j++)
        {
            DrawingUtilNG::hsv2rgb(Panel[i][j].pressure / 1000, 1., 1., red, green, blue);
            glColor3ub(red * 256, green * 256, blue * 256);
            if (i % 2 == 0)
                glVertex2f(i*20 , j*20);
            else
                glVertex2f(i*20+5, j*20+5);

            
        }
    }
    glEnd();
    glFlush();
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

    while (time < end_time)
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
        //Render(Panel);
        time += dt;
        //FsSleep(1);
    }
   
    return 0;
}
