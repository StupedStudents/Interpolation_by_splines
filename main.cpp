#include <QCoreApplication>
#include <QMatrix4x4>
#include <cmath>
#include <QDebug>
#include <fstream>

const int n=6;
double f(double x);

void sweep(float **matrix, float* b);

int main(int argc, char *argv[])
{
    double *points_x = new double[n];
    float **matrix = new float*[n];
    for(int i=0;i<n;i++){
        matrix[i]=new float[n];
    }
    double *points_y = new double[n];
    float *b = new float[n-2];
    double *steps = new double[n-1];
    double *p_c = new double[n];
    double *p_a = new double[n-1];
    double *p_d = new double[n-1];
    double *p_b = new double[n-1];

    if(n==6){
        points_x[0]=-2;
        points_x[1]=-1;
        points_x[2]=0;
        points_x[3]=1;
        points_x[4]=3;
        points_x[5]=4;
    }else{
        points_x[0]=-2;
        points_x[1]=-1.5;
        points_x[2]=-1;
        points_x[3]=-0.5;
        points_x[4]=0;
        points_x[5]=0.5;
        points_x[6]=1;
        points_x[7]=2;
        points_x[8]=3;
        points_x[9]=3.5;
        points_x[10]=4;
    }


    for(int i=0;i<n;i++){
        points_y[i]=f(points_x[i]);
    }
    for(int i = 0; i < n-1; i++){
        p_a[i] = points_y[i];
    }

    for(int i = 0; i < n-1; i++){
        steps[i]=points_x[i+1]-points_x[i];
    }
    matrix[0][0]=2*(steps[0]+steps[1]);
    matrix[0][1]=steps[1];
    matrix[n-3][n-3] = 2 * (steps[n-3] + steps[n - 2] );\
    matrix[n-3][n-2] = steps[n-3];
    for(int i = 1; i < n-2; i++)
    {
         matrix[i][i]=2*(steps[i]+steps[i+1]);
         matrix[i][i-1] = steps[i];
         matrix[i][i + 1] = steps[i + 1];
    }

    for(int i = 1; i < n-1; i++){
        b[i - 1] = 3 * ((points_y[i + 1] - points_y[i])/steps[i] - (points_y[i] - points_y[i - 1])/steps[i-1]);
    }

    sweep(matrix,b);

    p_c[0] = 0;
    p_c[n-1] = 0;
    for(int i = 1; i < n-1; i++){
        p_c[i] = b[i - 1];
    }
    for(int i = 0; i < n-1; i++){
        p_d[i] = (p_c[i + 1] - p_c[i])/ (3 * steps[i]);
    }
    for(int i = 0; i < n-1; i++){
        p_b[i] = (points_y[i + 1] - points_y[i]) / steps[i] - 1./3. * steps[i] * (p_c[i + 1] + 2 * p_c[i]);
    }

    std::ofstream STR("lab.txt");

    for(int i = 0; i < n-1; i++){
        float step = steps[i] / 100.;
        for(int j = 1; j < 101; j++){
            float tst = p_a[i] + p_b[i] * pow(step * j, 1) + p_c[i] * pow(step * j,2) + p_d[i] * pow(step * j, 3);
            STR  << (step * j + points_x[i])<< "\t" << tst;
            STR << "\n";
        }
        STR << "\n";
    }
    STR.close();
    return 0;
}

double inline f(double x){
    return (-2*pow(x,4) - 2*pow(x,3) - 2*pow(x,2) - 4*x - 2);
}

void sweep(float** a, float* b ){
    int N=n-2;
    int i;
    double znam;

    b[0]/=a[0][0];//Q1
    a[0][1]/=-a[0][0];//P1

    for(i=1;i < N-1;i++)
    {
    znam=-a[i][i]-a[i][i-1]*a[i-1][i]; //общий знаменатель для формул нахождения Pi, Qi
    a[i][i+1]/=znam; //Pi
    b[i]=(a[i][i-1]*b[i-1]-b[i])/znam; //Qi

    }
    //строка ниже для вычисления QN
    b[N-1]=(a[N-1][N-2]*b[N-2]-b[N-1])/(-a[N-1][N-1]-a[N-1][N-2]*a[N-2][N-1]);


    //обратный ход
    for(i=N-2;i > -1;i--)
    {
    b[i]+=b[i+1]*a[i][i+1];
    }


    return;
}
