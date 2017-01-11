#include <GL/glut.h>
#include <iostream>
#include <math.h>
#include <unistd.h>//WARNING: only used on linux OS, not working on windows.

GLsizei winWidth = 700,winHeight=700;//the size of the window

const static int total = 100;//the total times of the loop of the animating process

static int cycle = -5;//the current loop time,-5~-1 is used to get the four vertexs
const float lamda = 1.2;//the factor of the circle radius
static bool flag = 1;//used the choose between two modes
//the four vertexs of the quadralateral
static GLfloat
    PointA[2] ,
    PointB[2] ,
    PointC[2] ,
    PointD[2] ;

float max(float a,float b){
    return (a > b)? a : b ;
}

/*used to modify the number of cut-off points*/
int grad( int cycle){
    return 15-(int)sqrt(sqrt(cycle));
}

/*calculate the distance between two dots*/
float d(GLfloat p[],GLfloat q[]){
    return sqrt((p[0]-q[0])*(p[0]-q[0])+(p[1]-q[1])*(p[1]-q[1]));
}

/*plot a line with two points*/
void plotLine(GLfloat A[],GLfloat B[]){
    glBegin(GL_LINES);
    glVertex2fv(A);
    glVertex2fv(B);
    glEnd();
}
/*return x of the public point of two lines*/
GLfloat itsecX(GLfloat A[],GLfloat B[],GLfloat C[],GLfloat D[]){
    float a1 = A[0], a2 = A[1],
          b1 = B[0], b2 = B[1],
          c1 = C[0], c2 = C[1],
          d1 = D[0], d2 = D[1];

    float k1 = (b2-a2)/(b1-a1),
          k2 = (d2-c2)/(d1-c1);
    GLfloat X = (b2-d2+k2*d1-k1*b1)/(k2-k1);
    return X;
}
/*return y of the public point of two lines*/
GLfloat itsecY(GLfloat A[],GLfloat B[],GLfloat C[],GLfloat D[]){
    float a1 = A[0], a2 = A[1],
          b1 = B[0], b2 = B[1],
          c1 = C[0], c2 = C[1],
          d1 = D[0], d2 = D[1];

    float t1 = (b1-a1)/(b2-a2),
          t2 = (d1-c1)/(d2-c2);
    GLfloat Y = (b1-d1+t2*d2-t1*b2)/(t2-t1);
    return Y;
}
/*verify whether two points are the same one*/
bool ispoint(GLfloat A[],GLfloat B[]){
    bool
    a = (A[0]-B[0]>-0.0000001) && (A[0]-B[0]<0.0000001),
    b = (A[1]-B[1]>-0.0000001) && (A[1]-B[1]<0.0000001);
    return a&&b;
}
/*draw the delauny line in the giving sub-quadrilateral*/
void delauny(GLfloat A[],GLfloat B[],GLfloat C[],GLfloat D[]){
    if(ispoint(A,B)||ispoint(C,B)||ispoint(D,B)||ispoint(A,C)||ispoint(A,D)||ispoint(D,C))
        return;
    else{
        float
        AB = d(A,B),
        CB = d(B,C),
        CD = d(C,D),
        AD = d(D,A),
        BD = d(B,D);

        float
        sinA,cosA,sinC,cosC,sin_A_plus_C;

        cosA = (AB*AB+AD*AD-BD*BD)/(2*AB*AD);
        cosC = (CB*CB+CD*CD-BD*BD)/(2*CB*CD);
        sinA = sqrt(1-cosA*cosA);
        sinC = sqrt(1-cosC*cosC);
        sin_A_plus_C = sinA*cosC+cosA*sinC;

        if(100*sin_A_plus_C)
            plotLine(A,C);
        else
            plotLine(B,D);
    }
}


/*main display module*/
void myDisplay(void)
{
    glClear(GL_COLOR_BUFFER_BIT);

    if(cycle < 0){
        glColor3f(1.0f,0.0f,0.0f);
        std::cout << "cycle = " <<cycle << std::endl;
        if(cycle == -4){
            glBegin(GL_POINT);
            glPointSize(3.0f);
            glVertex2fv(PointA);
            glEnd();
        }
        else if(cycle == -3)
            plotLine(PointA,PointB);
        else if(cycle == -2){
            plotLine(PointA,PointB);
            plotLine(PointB,PointC);
        }
        else if(cycle == -1){
            plotLine(PointA,PointB);
            plotLine(PointB,PointC);
            plotLine(PointD,PointC);
            plotLine(PointD,PointA);
            cycle++;
        }
        glFlush();
        glutSwapBuffers();
        return;
    }
    if(cycle == 1){
        std::cout << "Choose the mode: 0 for pane ; 1 for ray." << std::endl;
//        std::cin >> flag;
        flag = 0;
    }
    GLfloat PointO[2],PointAB[2],PointBC[2],PointCD[2],PointDA[2];
    for(int i = 0;i < 2;i++){
        PointO[i] = 0.25*(PointA[i]+PointB[i]+PointC[i]+PointD[i]);
        PointAB[i] = 0.5*(PointA[i]+PointB[i]);
        PointBC[i] = 0.5*(PointC[i]+PointB[i]);
        PointCD[i] = 0.5*(PointC[i]+PointD[i]);
        PointDA[i] = 0.5*(PointA[i]+PointD[i]);
    }
    float R = lamda*max(max(d(PointO,PointA),d(PointO,PointB)),max(d(PointO,PointC),d(PointO,PointD)));
    float l,m;//used as factors
    GLfloat PointT[2];
    GLfloat PointX[2],PointY[2],PointP[2],PointQ[2];
    GLfloat PointX1[2],PointX2[2],PointY1[2],PointY2[2];
    GLfloat PointT1[2],PointT2[2],PointT3[2],PointT4[2];

    //draw the outline
    GLfloat r,g,b;
    if (cycle>=0&&cycle<=total/5){
        r=255;
        g=cycle*255/(total/5);
        b=0;
    }
    else if (cycle>total/5&&cycle<=2*total/5){
        r=255-(cycle-total/5)*255/(total/5);
        g=255;
        b=0;
    }
    else if(cycle>2*total/5&&cycle<=3*total/5){
        r=0;
        g=255;
        b=(cycle-2*total/5)*255/(total/5);
    }
    else if(cycle>3*total/5&&cycle<=4*total/5){
        r=0;
        g=255-(cycle-3*total/5)*255/(total/5);
        b=255;
    }
    else{
        r=(cycle-4*total/5)*255/(total/5);
        g=0;
        b=255;
    }
    glColor3f(r/255.0f,g/255.0f,b/255.0f);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    if(flag == 1){
        int n = 3;
        float nt = n;
        int no = (int)1.414*nt;
        int n_ab = (int)nt*sqrt(d(PointA,PointB))/d(PointO,PointAB);
        int n_bc = (int)nt*sqrt(d(PointB,PointC))/d(PointO,PointBC);
        int n_cd = (int)nt*sqrt(d(PointC,PointD))/d(PointO,PointCD);
        int n_da = (int)nt*sqrt(d(PointD,PointA))/d(PointO,PointDA);
        std::cout << "\r" << "cycle = " << cycle
              << "\tn=" << n
              << "\tn_ab=" << n_ab
              << "\tn_bc=" << n_bc
              << "\tn_cd=" << n_cd
              << "\tn_da=" << n_da
              << " d = "<<d(PointO,PointAB)
              << std::endl;
        for(int k = 0;k<no;k++){

            glBegin(GL_POLYGON);
            for(int i = 0;i<n_ab;i++){
                for(int j = 0;j < 2;j++)
                    PointT[j] = ((n_ab-i)*PointA[j]+i*PointB[j])/n_ab;
                l = d(PointT,PointO);
                m = (R-l)/total*cycle;
                for(int j = 0;j < 2;j++)
                    PointT[j]=((l+m)*PointT[j]-m*PointO[j])/l;
                for(int j = 0;j < 2;j++)
                    PointT[j]=(k*PointO[j]+(no-k)*PointT[j])/no;
                glVertex2fv(PointT);
            }

            for(int i = 0;i<n_bc;i++){
                for(int j = 0;j < 2;j++)
                    PointT[j] = ((n_bc-i)*PointB[j]+i*PointC[j])/n_bc;
                l = d(PointT,PointO);
                m = (R-l)/total*cycle;
                for(int j = 0;j < 2;j++)
                    PointT[j]=((l+m)*PointT[j]-m*PointO[j])/l;
                for(int j = 0;j < 2;j++)
                    PointT[j]=(k*PointO[j]+(no-k)*PointT[j])/no;
                glVertex2fv(PointT);
            }

            for(int i = 0;i<n_cd;i++){
                for(int j = 0;j < 2;j++)
                    PointT[j] = ((n_cd-i)*PointC[j]+i*PointD[j])/n_cd;
                l = d(PointT,PointO);
                m = (R-l)/total*cycle;
                for(int j = 0;j < 2;j++)
                    PointT[j]=((l+m)*PointT[j]-m*PointO[j])/l;
                for(int j = 0;j < 2;j++)
                    PointT[j]=(k*PointO[j]+(no-k)*PointT[j])/no;
                glVertex2fv(PointT);
            }

            for(int i = 0;i<n_da;i++){
                for(int j = 0;j < 2;j++)
                    PointT[j] = ((n_da-i)*PointD[j]+i*PointA[j])/n_da;
                l = d(PointT,PointO);
                m = (R-l)/total*cycle;
                for(int j = 0;j < 2;j++)
                    PointT[j]=((l+m)*PointT[j]-m*PointO[j])/l;
                for(int j = 0;j < 2;j++)
                    PointT[j]=(k*PointO[j]+(no-k)*PointT[j])/no;
                glVertex2fv(PointT);
            }
            glEnd();
        }

        //draw the sub-quadrilaterals
        for(int i = 0;i < n_ab+1;i++){
            for(int j = 0;j < 2;j++)
                PointX[j] = ((n_ab-i)*PointA[j]+i*PointB[j])/n_ab;

            l = d(PointX,PointO);
            m = (R-l)/total*cycle;
            for(int j = 0;j < 2;j++)
                PointX[j]=((l+m)*PointX[j]-m*PointO[j])/l;

            glBegin(GL_LINES);
            glVertex2fv(PointX);
            glVertex2fv(PointO);
            glEnd();
        }

        for(int i = 0;i < n_bc+1;i++){
            for(int j = 0;j < 2;j++)
                PointY[j] = ((n_bc-i)*PointB[j]+i*PointC[j])/n_bc;
            l = d(PointY,PointO);
            m = (R-l)/total*cycle;
            for(int j = 0;j < 2;j++)
                PointY[j]=((l+m)*PointY[j]-m*PointO[j])/l;

            glBegin(GL_LINES);
            glVertex2fv(PointY);
            glVertex2fv(PointO);
            glEnd();
        }

        for(int i = 0;i < n_cd+1;i++){
            for(int j = 0;j < 2;j++)
                PointP[j] = ((n_cd-i)*PointC[j]+i*PointD[j])/n_cd;
            l = d(PointP,PointO);
            m = (R-l)/total*cycle;
            for(int j = 0;j < 2;j++)
                PointP[j]=((l+m)*PointP[j]-m*PointO[j])/l;
            glBegin(GL_LINES);
            glVertex2fv(PointP);
            glVertex2fv(PointO);
            glEnd();

        }

        for(int i = 0;i < n_da+1;i++){
            for(int j = 0;j < 2;j++)
                PointQ[j] = ((n_da-i)*PointD[j]+i*PointA[j])/n_da;
            l = d(PointQ,PointO);
            m = (R-l)/total*cycle;
            for(int j = 0;j < 2;j++)
                PointQ[j]=((l+m)*PointQ[j]-m*PointO[j])/l;

            glBegin(GL_LINES);
            glVertex2fv(PointQ);
            glVertex2fv(PointO);
            glEnd();

        }

        //delauny part
        //AB
        for(int i = 0;i < n_ab;i++){

            for(int j = 0;j < 2;j++)
                PointX[j] = ((n_ab-i)*PointA[j]+i*PointB[j])/n_ab;
            l = d(PointX,PointO);
            m = (R-l)/total*cycle;
            for(int j = 0;j < 2;j++)
                PointX[j]=((l+m)*PointX[j]-m*PointO[j])/l;

            for(int j = 0;j < 2;j++)
                PointY[j] = ((n_ab-i-1)*PointA[j]+(i+1)*PointB[j])/n_ab;
            l = d(PointY,PointO);
            m = (R-l)/total*cycle;
            for(int j = 0;j < 2;j++)
                PointY[j]=((l+m)*PointY[j]-m*PointO[j])/l;

            for(int k = 1;k < no;k++){

                for(int j = 0;j < 2;j++)
                    PointX1[j]=(k*PointX[j]+(no-k)*PointO[j])/no;
                for(int j = 0;j < 2;j++)
                    PointX2[j]=((k+1)*PointX[j]+(no-k-1)*PointO[j])/no;
                for(int j = 0;j < 2;j++)
                    PointY1[j]=(k*PointY[j]+(no-k)*PointO[j])/no;
                for(int j = 0;j < 2;j++)
                    PointY2[j]=((k+1)*PointY[j]+(no-k-1)*PointO[j])/no;

                delauny(PointX1,PointX2,PointY2,PointY1);
            }

        }
        //BC
        for(int i = 0;i < n_bc;i++){

            for(int j = 0;j < 2;j++)
                PointX[j] = ((n_bc-i)*PointB[j]+i*PointC[j])/n_bc;
            l = d(PointX,PointO);
            m = (R-l)/total*cycle;
            for(int j = 0;j < 2;j++)
                PointX[j]=((l+m)*PointX[j]-m*PointO[j])/l;

            for(int j = 0;j < 2;j++)
                PointY[j] = ((n_bc-i-1)*PointB[j]+(i+1)*PointC[j])/n_bc;
            l = d(PointY,PointO);
            m = (R-l)/total*cycle;
            for(int j = 0;j < 2;j++)
                PointY[j]=((l+m)*PointY[j]-m*PointO[j])/l;

            for(int k = 1;k < no;k++){

                for(int j = 0;j < 2;j++)
                    PointX1[j]=(k*PointX[j]+(no-k)*PointO[j])/no;
                for(int j = 0;j < 2;j++)
                    PointX2[j]=((k+1)*PointX[j]+(no-k-1)*PointO[j])/no;
                for(int j = 0;j < 2;j++)
                    PointY1[j]=(k*PointY[j]+(no-k)*PointO[j])/no;
                for(int j = 0;j < 2;j++)
                    PointY2[j]=((k+1)*PointY[j]+(no-k-1)*PointO[j])/no;

                delauny(PointX1,PointX2,PointY2,PointY1);
            }
        }
        //CD
        for(int i = 0;i < n_cd;i++){

            for(int j = 0;j < 2;j++)
                PointX[j] = ((n_cd-i)*PointC[j]+i*PointD[j])/n_cd;
            l = d(PointX,PointO);
            m = (R-l)/total*cycle;
            for(int j = 0;j < 2;j++)
                PointX[j]=((l+m)*PointX[j]-m*PointO[j])/l;

            for(int j = 0;j < 2;j++)
                PointY[j] = ((n_cd-i-1)*PointC[j]+(i+1)*PointD[j])/n_cd;
            l = d(PointY,PointO);
            m = (R-l)/total*cycle;
            for(int j = 0;j < 2;j++)
                PointY[j]=((l+m)*PointY[j]-m*PointO[j])/l;

            for(int k = 1;k < no;k++){

                for(int j = 0;j < 2;j++)
                    PointX1[j]=(k*PointX[j]+(no-k)*PointO[j])/no;
                for(int j = 0;j < 2;j++)
                    PointX2[j]=((k+1)*PointX[j]+(no-k-1)*PointO[j])/no;
                for(int j = 0;j < 2;j++)
                    PointY1[j]=(k*PointY[j]+(no-k)*PointO[j])/no;
                for(int j = 0;j < 2;j++)
                    PointY2[j]=((k+1)*PointY[j]+(no-k-1)*PointO[j])/no;

                delauny(PointX1,PointX2,PointY2,PointY1);
            }

        }
        //DA
        for(int i = 0;i < n_da;i++){

            for(int j = 0;j < 2;j++)
                PointX[j] = ((n_da-i)*PointD[j]+i*PointA[j])/n_da;
            l = d(PointX,PointO);
            m = (R-l)/total*cycle;
            for(int j = 0;j < 2;j++)
                PointX[j]=((l+m)*PointX[j]-m*PointO[j])/l;

            for(int j = 0;j < 2;j++)
                PointY[j] = ((n_da-i-1)*PointD[j]+(i+1)*PointA[j])/n_da;
            l = d(PointY,PointO);
            m = (R-l)/total*cycle;
            for(int j = 0;j < 2;j++)
                PointY[j]=((l+m)*PointY[j]-m*PointO[j])/l;

            for(int k = 1;k < no;k++){

                for(int j = 0;j < 2;j++)
                    PointX1[j]=(k*PointX[j]+(no-k)*PointO[j])/no;
                for(int j = 0;j < 2;j++)
                    PointX2[j]=((k+1)*PointX[j]+(no-k-1)*PointO[j])/no;
                for(int j = 0;j < 2;j++)
                    PointY1[j]=(k*PointY[j]+(no-k)*PointO[j])/no;
                for(int j = 0;j < 2;j++)
                    PointY2[j]=((k+1)*PointY[j]+(no-k-1)*PointO[j])/no;

                delauny(PointX1,PointX2,PointY2,PointY1);
            }

        }

        int temp = 4*no;
        l = d(PointA,PointO);
        m = (R-l)/total*cycle;
        for(int j = 0;j < 2;j++)
            PointT1[j]=((l+m)*PointA[j]-m*PointO[j])/l;
        for(int j = 0;j < 2;j++)
            PointT1[j]=(PointT1[j]+(temp-1)*PointO[j])/(temp);

        l = d(PointB,PointO);
        m = (R-l)/total*cycle;
        for(int j = 0;j < 2;j++)
            PointT2[j]=((l+m)*PointB[j]-m*PointO[j])/l;
        for(int j = 0;j < 2;j++)
            PointT2[j]=(PointT2[j]+(temp-1)*PointO[j])/(temp);

        l = d(PointC,PointO);
        m = (R-l)/total*cycle;
        for(int j = 0;j < 2;j++)
            PointT3[j]=((l+m)*PointC[j]-m*PointO[j])/l;
        for(int j = 0;j < 2;j++)
            PointT3[j]=(PointT3[j]+(temp-1)*PointO[j])/(temp);

        l = d(PointD,PointO);
        m = (R-l)/total*cycle;
        for(int j = 0;j < 2;j++)
            PointT4[j]=((l+m)*PointD[j]-m*PointO[j])/l;
        for(int j = 0;j < 2;j++)
            PointT4[j]=(PointT4[j]+(temp-1)*PointO[j])/(temp);

        glColor3f(0.0f,0.0f,0.0f);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glBegin(GL_POLYGON);
        glVertex2fv(PointT1);
        glVertex2fv(PointT2);
        glVertex2fv(PointT3);
        glVertex2fv(PointT4);
        glEnd();

        glColor3f(r/255.0f,g/255.0f,b/255.0f);
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glBegin(GL_POLYGON);
        glVertex2fv(PointT1);
        glVertex2fv(PointT2);
        glVertex2fv(PointT3);
        glVertex2fv(PointT4);
        glEnd();

        delauny(PointT1,PointT2,PointT3,PointT4);
    }

    else if(flag == 0){
        int n = grad(cycle);
        std::cout << "\r" << "cycle = " << cycle
              << "\t\tn=" << n << std::endl;
        //draw the outline
        glBegin(GL_POLYGON);
        GLfloat PointT[2];

        for(int i = 0;i<n+1;i++){
            for(int j = 0;j < 2;j++)
                PointT[j] = ((n-i)*PointA[j]+i*PointB[j])/n;
            l = d(PointT,PointO);
            m = (R-l)/total*cycle;
            for(int j = 0;j < 2;j++)
                PointT[j]=((l+m)*PointT[j]-m*PointO[j])/l;
            glVertex2fv(PointT);
        }

        for(int i = 0;i<n+1;i++){
            for(int j = 0;j < 2;j++)
                PointT[j] = ((n-i)*PointB[j]+i*PointC[j])/n;
            l = d(PointT,PointO);
            m = (R-l)/total*cycle;
            for(int j = 0;j < 2;j++)
                PointT[j]=((l+m)*PointT[j]-m*PointO[j])/l;
            glVertex2fv(PointT);
        }

        for(int i = 0;i<n+1;i++){
            for(int j = 0;j < 2;j++)
                PointT[j] = ((n-i)*PointC[j]+i*PointD[j])/n;
            l = d(PointT,PointO);
            m = (R-l)/total*cycle;
            for(int j = 0;j < 2;j++)
                PointT[j]=((l+m)*PointT[j]-m*PointO[j])/l;
            glVertex2fv(PointT);
        }

        for(int i = 0;i<n+1;i++){
            for(int j = 0;j < 2;j++)
                PointT[j] = ((n-i)*PointD[j]+i*PointA[j])/n;
            l = d(PointT,PointO);
            m = (R-l)/total*cycle;
            for(int j = 0;j < 2;j++)
                PointT[j]=((l+m)*PointT[j]-m*PointO[j])/l;
            glVertex2fv(PointT);
        }
        glEnd();


        //draw the sub-quadrilaterals
        GLfloat PointX[2],PointY[2],PointP[2],PointQ[2],
                PointX_[2],PointY_[2],PointP_[2],PointQ_[2];
        for(int i = 0;i < n+1;i++){
            for(int j = 0;j < 2;j++){
                PointX[j] = ((n-i)*PointA[j]+i*PointB[j])/n;
                PointY[j] = ((n-i)*PointD[j]+i*PointC[j])/n;
                PointP[j] = ((n-i)*PointA[j]+i*PointD[j])/n;
                PointQ[j] = ((n-i)*PointB[j]+i*PointC[j])/n;
            }

            l = d(PointX,PointO);
            m = (R-l)/total*cycle;
            for(int j = 0;j < 2;j++)
                PointX[j]=((l+m)*PointX[j]-m*PointO[j])/l;

            l = d(PointY,PointO);
            m = (R-l)/total*cycle;
            for(int j = 0;j < 2;j++)
                PointY[j]=((l+m)*PointY[j]-m*PointO[j])/l;

            l = d(PointP,PointO);
            m = (R-l)/total*cycle;
            for(int j = 0;j < 2;j++)
                PointP[j]=((l+m)*PointP[j]-m*PointO[j])/l;

            l = d(PointQ,PointO);
            m = (R-l)/total*cycle;
            for(int j = 0;j < 2;j++)
                PointQ[j]=((l+m)*PointQ[j]-m*PointO[j])/l;

            glBegin(GL_LINES);
            glVertex2fv(PointX);
            glVertex2fv(PointY);
            glEnd();

            glBegin(GL_LINES);
            glVertex2fv(PointP);
            glVertex2fv(PointQ);
            glEnd();


            if(i < n/2){
                for(int j = 0;j < 2;j++){
                PointX_[j] = ((n-i)*PointB[j]+i*PointA[j])/n;
                PointY_[j] = ((n-i)*PointC[j]+i*PointD[j])/n;
                PointP_[j] = ((n-i)*PointD[j]+i*PointA[j])/n;
                PointQ_[j] = ((n-i)*PointC[j]+i*PointB[j])/n;
                }

                l = d(PointX_,PointO);
                m = (R-l)/total*cycle;
                for(int j = 0;j < 2;j++)
                    PointX_[j]=((l+m)*PointX_[j]-m*PointO[j])/l;

                l = d(PointY_,PointO);
                m = (R-l)/total*cycle;
                for(int j = 0;j < 2;j++)
                    PointY_[j]=((l+m)*PointY_[j]-m*PointO[j])/l;

                l = d(PointP_,PointO);
                m = (R-l)/total*cycle;
                for(int j = 0;j < 2;j++)
                    PointP_[j]=((l+m)*PointP_[j]-m*PointO[j])/l;

                l = d(PointQ_,PointO);
                m = (R-l)/total*cycle;
                for(int j = 0;j < 2;j++)
                    PointQ_[j]=((l+m)*PointQ_[j]-m*PointO[j])/l;

                glBegin(GL_LINES);
                glVertex2fv(PointX);
                glVertex2fv(PointX_);
                glEnd();

                glBegin(GL_LINES);
                glVertex2fv(PointY);
                glVertex2fv(PointY_);
                glEnd();

                glBegin(GL_LINES);
                glVertex2fv(PointP);
                glVertex2fv(PointP_);
                glEnd();

                glBegin(GL_LINES);
                glVertex2fv(PointQ);
                glVertex2fv(PointQ_);
                glEnd();

            }

        }
        //draw the delauny triangle seperation on each sub-quadrilateral
        GLfloat PointX1[2],PointY1[2],PointP1[2],PointQ1[2],
                PointX2[2],PointY2[2],PointP2[2],PointQ2[2],
                PointX3[2],PointY3[2],PointP3[2],PointQ3[2],
                PointX4[2],PointY4[2],PointP4[2],PointQ4[2],
                PointX11[2],PointY11[2],PointP11[2],PointQ11[2],
                PointX22[2],PointY22[2],PointP22[2],PointQ22[2],
                PointT1[2],PointT2[2],PointT3[2],PointT4[2];
        for(int i = 0;i < n;i++){
            for(int j = 0;j < 2;j++){
                PointX1[j] = ((n-i)*PointB[j]+i*PointA[j])/n;
                PointX2[j] = ((n-(i+1))*PointB[j]+(i+1)*PointA[j])/n;
                PointY1[j] = ((n-i)*PointC[j]+i*PointD[j])/n;
                PointY2[j] = ((n-(i+1))*PointC[j]+(i+1)*PointD[j])/n;
            }

            l = d(PointX1,PointO);
            m = (R-l)/total*cycle;
            for(int j = 0;j < 2;j++)
                PointX1[j]=((l+m)*PointX1[j]-m*PointO[j])/l;

            l = d(PointX2,PointO);
            m = (R-l)/total*cycle;
            for(int j = 0;j < 2;j++)
                PointX2[j]=((l+m)*PointX2[j]-m*PointO[j])/l;

            l = d(PointY1,PointO);
            m = (R-l)/total*cycle;
            for(int j = 0;j < 2;j++)
                PointY1[j]=((l+m)*PointY1[j]-m*PointO[j])/l;

            l = d(PointY2,PointO);
            m = (R-l)/total*cycle;
            for(int j = 0;j < 2;j++)
                PointY2[j]=((l+m)*PointY2[j]-m*PointO[j])/l;

            if(i<n/2){
                for(int j = 0;j < 2;j++){
                    PointX11[j] = ((n-i)*PointA[j]+i*PointB[j])/n;
                    PointX22[j] = ((n-(i+1))*PointA[j]+(i+1)*PointB[j])/n;
                    PointY11[j] = ((n-i)*PointD[j]+i*PointC[j])/n;
                    PointY22[j] = ((n-(i+1))*PointD[j]+(i+1)*PointC[j])/n;
                }

                l = d(PointX11,PointO);
                m = (R-l)/total*cycle;
                for(int j = 0;j < 2;j++)
                    PointX11[j]=((l+m)*PointX11[j]-m*PointO[j])/l;

                l = d(PointX22,PointO);
                m = (R-l)/total*cycle;
                for(int j = 0;j < 2;j++)
                    PointX22[j]=((l+m)*PointX22[j]-m*PointO[j])/l;

                l = d(PointY11,PointO);
                m = (R-l)/total*cycle;
                for(int j = 0;j < 2;j++)
                    PointY11[j]=((l+m)*PointY11[j]-m*PointO[j])/l;

                l = d(PointY22,PointO);
                m = (R-l)/total*cycle;
                for(int j = 0;j < 2;j++)
                    PointY22[j]=((l+m)*PointY22[j]-m*PointO[j])/l;

                for(int ii = i; ii < n-2-i ;ii++){
                    for(int j = 0;j < 2;j++){
                        PointX4[j] = ((n-(ii+1))*PointB[j]+(ii+1)*PointA[j])/n;
                        PointY4[j] = ((n-(ii+1))*PointC[j]+(ii+1)*PointD[j])/n;
                        PointX3[j] = ((n-(ii+2))*PointB[j]+(ii+2)*PointA[j])/n;
                        PointY3[j] = ((n-(ii+2))*PointC[j]+(ii+2)*PointD[j])/n;
                    }
                    l = d(PointX3,PointO);
                    m = (R-l)/total*cycle;
                    for(int j = 0;j < 2;j++)
                        PointX3[j]=((l+m)*PointX3[j]-m*PointO[j])/l;

                    l = d(PointY3,PointO);
                    m = (R-l)/total*cycle;
                    for(int j = 0;j < 2;j++)
                        PointY3[j]=((l+m)*PointY3[j]-m*PointO[j])/l;

                    l = d(PointX4,PointO);
                    m = (R-l)/total*cycle;
                    for(int j = 0;j < 2;j++)
                        PointX4[j]=((l+m)*PointX4[j]-m*PointO[j])/l;

                    l = d(PointY4,PointO);
                    m = (R-l)/total*cycle;
                    for(int j = 0;j < 2;j++)
                        PointY4[j]=((l+m)*PointY4[j]-m*PointO[j])/l;


                    PointT1[0] = itsecX(PointX1,PointX11,PointX4,PointY4);
                    PointT1[1] = itsecY(PointX1,PointX11,PointX4,PointY4);
                    PointT2[0] = itsecX(PointX1,PointX11,PointX3,PointY3);
                    PointT2[1] = itsecY(PointX1,PointX11,PointX3,PointY3);
                    PointT3[0] = itsecX(PointX2,PointX22,PointX4,PointY4);
                    PointT3[1] = itsecY(PointX2,PointX22,PointX4,PointY4);
                    PointT4[0] = itsecX(PointX2,PointX22,PointX3,PointY3);
                    PointT4[1] = itsecY(PointX2,PointX22,PointX3,PointY3);
                    delauny(PointT1,PointT2,PointT4,PointT3);

                    PointT1[0] = itsecX(PointY1,PointY11,PointX4,PointY4);
                    PointT1[1] = itsecY(PointY1,PointY11,PointX4,PointY4);
                    PointT2[0] = itsecX(PointY1,PointY11,PointX3,PointY3);
                    PointT2[1] = itsecY(PointY1,PointY11,PointX3,PointY3);
                    PointT3[0] = itsecX(PointY2,PointY22,PointX4,PointY4);
                    PointT3[1] = itsecY(PointY2,PointY22,PointX4,PointY4);
                    PointT4[0] = itsecX(PointY2,PointY22,PointX3,PointY3);
                    PointT4[1] = itsecY(PointY2,PointY22,PointX3,PointY3);
                    delauny(PointT1,PointT2,PointT4,PointT3);
                }


            }


            for(int j = 0;j < n;j++){
                for(int k = 0;k < 2;k++){
                PointP1[k] = ((n-j)*PointB[k]+j*PointC[k])/n;
                PointP2[k] = ((n-(j+1))*PointB[k]+(j+1)*PointC[k])/n;
                PointQ1[k] = ((n-j)*PointA[k]+j*PointD[k])/n;
                PointQ2[k] = ((n-(j+1))*PointA[k]+(j+1)*PointD[k])/n;
                }

                l = d(PointP1,PointO);
                m = (R-l)/total*cycle;
                for(int j = 0;j < 2;j++)
                    PointP1[j]=((l+m)*PointP1[j]-m*PointO[j])/l;

                l = d(PointP2,PointO);
                m = (R-l)/total*cycle;
                for(int j = 0;j < 2;j++)
                    PointP2[j]=((l+m)*PointP2[j]-m*PointO[j])/l;

                l = d(PointQ1,PointO);
                m = (R-l)/total*cycle;
                for(int j = 0;j < 2;j++)
                    PointQ1[j]=((l+m)*PointQ1[j]-m*PointO[j])/l;

                l = d(PointQ2,PointO);
                m = (R-l)/total*cycle;
                for(int j = 0;j < 2;j++)
                    PointQ2[j]=((l+m)*PointQ2[j]-m*PointO[j])/l;

                PointT1[0] = itsecX(PointX1,PointY1,PointP1,PointQ1);
                PointT1[1] = itsecY(PointX1,PointY1,PointP1,PointQ1);
                PointT2[0] = itsecX(PointX1,PointY1,PointP2,PointQ2);
                PointT2[1] = itsecY(PointX1,PointY1,PointP2,PointQ2);
                PointT3[0] = itsecX(PointX2,PointY2,PointP1,PointQ1);
                PointT3[1] = itsecY(PointX2,PointY2,PointP1,PointQ1);
                PointT4[0] = itsecX(PointX2,PointY2,PointP2,PointQ2);
                PointT4[1] = itsecY(PointX2,PointY2,PointP2,PointQ2);
                delauny(PointT1,PointT2,PointT4,PointT3);

                if(i<n/2){
                    for(int k = 0;k < 2;k++){
                        PointP11[k] = ((n-j)*PointC[k]+j*PointB[k])/n;
                        PointP22[k] = ((n-(j+1))*PointC[k]+(j+1)*PointB[k])/n;
                        PointQ11[k] = ((n-j)*PointD[k]+j*PointA[k])/n;
                        PointQ22[k] = ((n-(j+1))*PointD[k]+(j+1)*PointA[k])/n;
                    }

                    l = d(PointP11,PointO);
                    m = (R-l)/total*cycle;
                    for(int j = 0;j < 2;j++)
                        PointP11[j]=((l+m)*PointP11[j]-m*PointO[j])/l;

                    l = d(PointP22,PointO);
                    m = (R-l)/total*cycle;
                    for(int j = 0;j < 2;j++)
                        PointP22[j]=((l+m)*PointP22[j]-m*PointO[j])/l;

                    l = d(PointQ11,PointO);
                    m = (R-l)/total*cycle;
                    for(int j = 0;j < 2;j++)
                        PointQ11[j]=((l+m)*PointQ11[j]-m*PointO[j])/l;

                    l = d(PointQ22,PointO);
                    m = (R-l)/total*cycle;
                    for(int j = 0;j < 2;j++)
                        PointQ22[j]=((l+m)*PointQ22[j]-m*PointO[j])/l;

                    for(int ii = j; ii < n-2-j;ii++){
                        for(int j = 0;j < 2;j++){
                            PointP4[j] = ((n-(ii+1))*PointB[j]+(ii+1)*PointC[j])/n;
                            PointQ4[j] = ((n-(ii+1))*PointA[j]+(ii+1)*PointD[j])/n;
                            PointP3[j] = ((n-(ii+2))*PointB[j]+(ii+2)*PointC[j])/n;
                            PointQ3[j] = ((n-(ii+2))*PointA[j]+(ii+2)*PointD[j])/n;
                        }
                        l = d(PointP3,PointO);
                        m = (R-l)/total*cycle;
                        for(int j = 0;j < 2;j++)
                            PointP3[j]=((l+m)*PointP3[j]-m*PointO[j])/l;

                        l = d(PointQ3,PointO);
                        m = (R-l)/total*cycle;
                        for(int j = 0;j < 2;j++)
                            PointQ3[j]=((l+m)*PointQ3[j]-m*PointO[j])/l;

                        l = d(PointP4,PointO);
                        m = (R-l)/total*cycle;
                        for(int j = 0;j < 2;j++)
                            PointP4[j]=((l+m)*PointP4[j]-m*PointO[j])/l;

                        l = d(PointQ4,PointO);
                        m = (R-l)/total*cycle;
                        for(int j = 0;j < 2;j++)
                            PointQ4[j]=((l+m)*PointQ4[j]-m*PointO[j])/l;

                        PointT1[0] = itsecX(PointP1,PointP11,PointP4,PointQ4);
                        PointT1[1] = itsecY(PointP1,PointP11,PointP4,PointQ4);
                        PointT2[0] = itsecX(PointP1,PointP11,PointP3,PointQ3);
                        PointT2[1] = itsecY(PointP1,PointP11,PointP3,PointQ3);
                        PointT3[0] = itsecX(PointP2,PointP22,PointP4,PointQ4);
                        PointT3[1] = itsecY(PointP2,PointP22,PointP4,PointQ4);
                        PointT4[0] = itsecX(PointP2,PointP22,PointP3,PointQ3);
                        PointT4[1] = itsecY(PointP2,PointP22,PointP3,PointQ3);
                        delauny(PointT1,PointT2,PointT4,PointT3);

                        PointT1[0] = itsecX(PointQ1,PointQ11,PointP4,PointQ4);
                        PointT1[1] = itsecY(PointQ1,PointQ11,PointP4,PointQ4);
                        PointT2[0] = itsecX(PointQ1,PointQ11,PointP3,PointQ3);
                        PointT2[1] = itsecY(PointQ1,PointQ11,PointP3,PointQ3);
                        PointT3[0] = itsecX(PointQ2,PointQ22,PointP4,PointQ4);
                        PointT3[1] = itsecY(PointQ2,PointQ22,PointP4,PointQ4);
                        PointT4[0] = itsecX(PointQ2,PointQ22,PointP3,PointQ3);
                        PointT4[1] = itsecY(PointQ2,PointQ22,PointP3,PointQ3);
                        delauny(PointT1,PointT2,PointT4,PointT3);
                    }
                }

            }

        }

    }
    glFlush();

    glutSwapBuffers();
}

// use the mouse click to get the four vertexes
void myMouse(int btn, int state, int x, int y)
{
    if(btn == GLUT_LEFT_BUTTON && state == GLUT_DOWN && cycle==-5){
        PointA[0] = (GLfloat)2*x/winWidth-1.0;PointA[1] = -((GLfloat)2*y/winHeight-1.0);
        std::cout << "the position of A is " << PointA[0] << " " << PointA[1] << std::endl;
        cycle++;
        myDisplay();
    }
    else if(btn == GLUT_LEFT_BUTTON && state == GLUT_DOWN && cycle==-4){
        PointB[0] = (GLfloat)2*x/winWidth-1.0;PointB[1] = -((GLfloat)2*y/winHeight-1.0);
        std::cout << "the position of B is " << PointB[0] << " " << PointB[1] << std::endl;
        cycle++;
        myDisplay();
    }
    else if(btn == GLUT_LEFT_BUTTON && state == GLUT_DOWN && cycle==-3){
        PointC[0] = (GLfloat)2*x/winWidth-1.0;PointC[1] = -((GLfloat)2*y/winHeight-1.0);
        std::cout << "the position of C is " << PointC[0] << " " << PointC[1] << std::endl;
        cycle++;
        myDisplay();
    }
    else if(btn == GLUT_LEFT_BUTTON && state == GLUT_DOWN && cycle==-2){
        PointD[0] = (GLfloat)2*x/winWidth-1.0;PointD[1] = -((GLfloat)2*y/winHeight-1.0);
        std::cout << "the position of D is " << PointD[0] << " " << PointD[1] << std::endl;
        cycle++;
        myDisplay();
    }
    else if(btn == GLUT_RIGHT_BUTTON && state == GLUT_DOWN){
        cycle = -5;
        myDisplay();
    }
    else if(btn == GLUT_MIDDLE_BUTTON && state == GLUT_DOWN){
        exit(0);
    }
    glFlush();
    return;
}

void myIdle(void)
{
     if(cycle<0)
        return;
     cycle++;
     if( cycle >= total ){
         cycle = -5;
         sleep(5);
     }
     myDisplay();
}

int main(int argc, char *argv[]){
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
    glOrtho(1.0,-1.0,-1.0,1.0,-1.0,1.0);
    glViewport(-1.0,-1.0,winWidth,winHeight);
    glutInitWindowPosition(200, 0);
    glutInitWindowSize(winWidth,winHeight);
    glutCreateWindow("Quadrilateral->Circle");
    glutMouseFunc(myMouse);
    glutDisplayFunc(&myDisplay);
    glutIdleFunc(&myIdle);
    glutMainLoop();
    system("pause");
    return 0;
}


