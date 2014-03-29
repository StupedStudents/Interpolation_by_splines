#include <QCoreApplication>
#include <QMatrix4x4>
#include <cmath>
#include <QDebug>

double f(double x);

int main(int argc, char *argv[])
{
    double *points_x = new double[6];
    points_x[0]=-1;
    points_x[1]=0;
    points_x[2]=1;
    points_x[3]=2;
    points_x[4]=3;
    points_x[5]=4;
    double *points_y = new double[6];
//    for(int i=0;i<6;i++){
//        points_y[i]=f(points_x[i]);
//    }

    points_y[0]=14;
    points_y[1]=3;
    points_y[2]=-2;
    points_y[3]=-7;
    points_y[4]=6;
    points_y[5]=79;

    QMatrix4x4 *matrix = new QMatrix4x4();
    double *b = new double[4];
    double *steps = new double[5];


    for(int i=0;i<5;i++){
        steps[i]=points_x[i+1]-points_x[i];
    }
    matrix->setRow(0,QVector4D(2*(steps[0]+steps[1]),steps[1],0,0));
    matrix->setRow(1,QVector4D(steps[1],2*(steps[1]+steps[2]),steps[2],0));
    matrix->setRow(2,QVector4D(0,steps[2],2*(steps[2]+steps[3]),steps[4]));
    matrix->setRow(3,QVector4D(0,0,steps[3],2*(steps[3]+steps[4])));
    for(int i = 1; i < 5; i++){
        b[i - 1] = 3 * ((points_y[i + 1] - points_y[i])/steps[i] - (points_y[i] - points_y[i - 1])/steps[i-1]);
    }


    return 0;
}

double inline f(double x){
    return (-2*pow(x,4) - 2*pow(x,3) - 2*pow(x,2) - 4*x - 2);
}

void sweep(QMatrix4x4 *matrix, double* b){
    double zam;

    b[0]/=a[0][0];//Q1
    a[0][1]/=-a[0][0];//P1

    for(int i=1;i < -1;i++){
        zam=-a[i][i]-a[i][i-1]*a[i-1][i]; //общий знаменатель для формул нахождения Pi, Qi
        a[i][i+1]/=zam; //Pi
        b[i]=(a[i][i-1]*b[i-1]-b[i])/zam; //Qi
    }
    //строка ниже для вычисления Q
    b[-1]=(a[-1][-2]*b[-2]-b[-1])/(-a[-1][-1]-a[-1][-2]*a[-2][-1]);

    //обратный ход
    for(int i=-2;i > -1;i--)
    {
    b[i]+=b[i+1]*a[i][i+1];
    }

    return;
}
