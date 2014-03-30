#include <QCoreApplication>
#include <QMatrix4x4>
#include <cmath>
#include <QDebug>
#include <fstream>

double f(double x);

void sweep(QMatrix4x4 *matrix, double* b);

int main(int argc, char *argv[])
{
    double *points_x = new double[6];
    points_x[0]=-2;
    points_x[1]=-1;
    points_x[2]=0;
    points_x[3]=1;
    points_x[4]=3;
    points_x[5]=4;
    double *points_y = new double[6];
    for(int i=0;i<6;i++){
        points_y[i]=f(points_x[i]);
    }

    //points_y[0]=14;
    //points_y[1]=3;
    //points_y[2]=-2;
    //points_y[3]=-7;
    //points_y[4]=6;
    //points_y[5]=79;

    QMatrix4x4 *matrix = new QMatrix4x4();
    double *b = new double[4];
    double *steps = new double[5];
    double *p_c = new double[6];
    double *p_a = new double[5];
    double *p_d = new double[5];
    double *p_b = new double[5];

    for(int i = 0; i < 5; i++)
    {
        p_a[i] = points_y[i];
    }

    for(int i = 0; i < 5; i++){
        steps[i]=points_x[i+1]-points_x[i];
    }
    matrix->setRow(0,QVector4D(2*(steps[0]+steps[1]),steps[1],0,0));
    matrix->setRow(1,QVector4D(steps[1],2*(steps[1]+steps[2]),steps[2],0));
    matrix->setRow(2,QVector4D(0,steps[2],2*(steps[2]+steps[3]),steps[4]));
    matrix->setRow(3,QVector4D(0,0,steps[3],2*(steps[3]+steps[4])));
    for(int i = 1; i < 5; i++){
        b[i - 1] = 3 * ((points_y[i + 1] - points_y[i])/steps[i] - (points_y[i] - points_y[i - 1])/steps[i-1]);
    }

    sweep(matrix,b);
    p_c[0] = 0;
    p_c[5] = 0;
    for(int i = 1; i < 5; i++)
    {
        p_c[i] = b[i - 1];
    }

    for(int i = 0; i < 5; i++)
    {
        p_d[i] = (p_c[i + 1] - p_c[i])/ (3 * steps[i]);
    }

    for(int i = 0; i < 5; i++)
    {
        p_b[i] = (points_y[i + 1] - points_y[i]) / steps[i] - 1./3. * steps[i] * (p_c[i + 1] + 2 * p_c[i]);
    }

    std::ofstream STR("lab.txt");

    for(int i = 0; i < 5; i++)
    {
        float step = steps[i] / 100.;
        for(int j = 1; j < 101; j++)
        {
            //float tst = p_a[i] * pow(step * j + points_x[i] , 3) + p_b[i] * pow(step * j + points_x[i], 2) + p_c[i] * (step * j + points_x[i]) + p_d[i] ;
            //float tst = p_d[i] * pow( step * j - points_x[i], 3) / 6. + p_c[i] * pow(step * j - points_x[i], 2) / 2.
              //      + p_b[i] * (step * j - points_x[i] ) + p_a[i] ;
            float tst = p_a[i] + p_b[i] * pow(step * j, 1) + p_c[i] * pow(step * j,2) + p_d[i] * pow(step * j, 3);
            STR << tst << "\t" << (step * j + points_x[i]);

            STR << "\n";
        }

        STR << "\n";
    }

    for(int i = 0; i < 6; i++)
    {
        qDebug() << points_y[i];
    }
    return 0;
}

double inline f(double x){
    return (-2*pow(x,4) - 2*pow(x,3) - 2*pow(x,2) - 4*x - 2);
}

void sweep(QMatrix4x4 *matrix, double* b){
    double zam;
    float *a = matrix->data();

    b[0] /= a[0];
    a[1] /= -a[0];

    for(int i = 1; i < 3; i++){
        zam = -a[i * 4 + i]-a[i * 4 + i - 1]*a[(i - 1) * 4 + i];
        a[i * 4 + i + 1] /= zam;
        b[i] = (a[i * 4 + i - 1] * b[i-1] - b[i]) / zam;
    }

    b[3] = (a[14] * b[2] - b[3]) / (-a[15] - a[14] * a[11]);

    for(int i = 2; i > -1; i--)
    {
        b[i] += b[i+1] * a[i * 4 + i + 1];
    }

    return;
}
