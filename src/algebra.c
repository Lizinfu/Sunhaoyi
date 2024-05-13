#include "algebra.h"
#include <stdio.h>

Matrix create_matrix(int row, int col)
{
    Matrix m;
    m.rows = row;
    m.cols = col;
    return m;
}

Matrix add_matrix(Matrix a, Matrix b)
{
    int i,j;
    if(a.rows==b.rows&&a.cols==b.cols){
        Matrix m = create_matrix(a.rows, a.cols);
        for (i = 0; i < a.rows; i++){
            for (j = 0; j < a.cols; j++){
    		    m.data[i][j]=a.data[i][j]+b.data[i][j];
    	    }
        }
        return m;
    }else{
        printf("Error: Matrix a and b must have the same rows and cols.");
    }
    return create_matrix(0, 0);
}

Matrix sub_matrix(Matrix a, Matrix b)
{
    int i,j;
    if(a.rows==b.rows&&a.cols==b.cols){
        Matrix m = create_matrix(a.rows, a.cols);
        for (i = 0; i < a.rows; i++){
            for (j = 0; j < a.cols; j++){
    		    m.data[i][j]=a.data[i][j]-b.data[i][j];
    	    }
        }
        return m;
    }else{
        printf("Error: Matrix a and b must have the same rows and cols.");
    }
    return create_matrix(0, 0);
}

Matrix mul_matrix(Matrix a, Matrix b)
{
    int i,j,r;
    if(a.cols==b.rows){
        Matrix m = create_matrix(a.rows, b.cols);
        for (i = 0; i < a.rows; i++){
            for(j = 0;j < b.cols; j++){
                m.data[i][j]=0;
            }
        }
        for (i = 0; i < a.rows; i++){
            for(j = 0;j < b.cols; j++){
                for(r = 0;r < a.cols;r++){
                    m.data[i][j]=m.data[i][j]+a.data[i][r]*b.data[r][j];
                }
            }
        }
        return m;
    }else{
        printf("Error: The number of cols of matrix a must be equal to the number of rows of matrix b.");
    }
    return create_matrix(0, 0);
}

Matrix scale_matrix(Matrix a, double k)
{
    int i,j;
    Matrix m = create_matrix(a.rows, a.cols);
    for (i = 0; i < a.rows; i++){
        for (j = 0; j < a.cols; j++){
    		m.data[i][j]=a.data[i][j]*k;
    	}
    }
    return m;
    return create_matrix(0, 0);
}

Matrix transpose_matrix(Matrix a)
{
    int i,j;
    Matrix m = create_matrix(a.cols, a.rows);
    for (i = 0; i < a.cols; i++){
        for (j = 0; j < a.rows; j++){
    		m.data[i][j]=a.data[j][i];
    	}
    }
    return m;
    return create_matrix(0, 0);
}

double det_matrix(Matrix a)
{
    double det=0;
    int i,r,p;
    if(a.cols==a.rows){
        if(a.cols==1){
            return a.data[0][0];
        }else if(a.cols>1){
            for(i=0;i<a.cols;i++){
                Matrix m=create_matrix(a.rows-1,a.cols-1);
                for(r=0;r<a.rows-1;r++){
                    for(p=0;p<a.cols;p++){
                        if(p<i){
                            m.data[r][p]=a.data[r+1][p];
                        }else if (p>=i){
                            m.data[r][p]=a.data[r+1][p+1];
                        }
                        
                    }
                }
                det=det+((i+1)%2==1?1:(-1))*det_matrix(m)*a.data[0][i];
            }
            return det;
        }
    }else{
        printf("Error: The matrix must be a square matrix.");
    }
    return 0;
}

Matrix inv_matrix(Matrix a)
{
   if(a.rows==a.cols&&det_matrix(a)!=0){
    	if(a.rows>1){
        int i,j,r,p;
        double det;
        Matrix m=create_matrix(a.rows,a.cols);
        for(i=0;i<a.rows;i++){
            for(j=0;j<a.cols;j++){
                Matrix n=create_matrix(a.rows-1,a.cols-1);
                for(r=0;r<a.rows-1;r++){
                    for(p=0;p<a.cols-1;p++){
                        if(r<i){
                            if(p<j){
                                n.data[r][p]=a.data[r][p];
                            }else if(p>=j){
                                n.data[r][p]=a.data[r][j+1];
                            }
                        }else if (r>=i){
                            if(p<j){
                                n.data[r][p]=a.data[r+1][p];
                            }else if(p>=j){
                                n.data[r][p]=a.data[r+1][p+1];
                            }
                        }
                        
                    }
                }
                m.data[i][j]=((i+j)%2==0?1:-1)*det_matrix(n);
            }
        }
        if(m.data[0]!=NULL){
            return m;
        }else{
            printf("Error: The matrix must be a square matrix.");
        }
        return create_matrix(0, 0);
    }else{
    	return a;
	}
    }else if(a.rows==a.cols&&det_matrix(a)==0){
        printf("Error: The matrix is singular.");
	return create_matrix(0, 0);
    }
}

int rank_matrix(Matrix a)
{
    int rank;
    int i,j,r;
    double t[a.cols];
    double k;
    rank=(a.rows>=a.cols?a.cols:a.rows);
    for(i=0;i<a.rows;i++){
        if(a.data[i][i]!=0){
        	for(r=i+1;r<a.rows;r++){
				k=a.data[r][i]/a.data[i][i];
        		for(j=i;j<a.cols;j++){
        			a.data[r][j]=a.data[r][j]-a.data[i][j]*k;
				}
			}
		}else{
			for(r=i+1;r<a.rows;r++){
				if(a.data[r][i]!=0){
					goto next;
				}
			}
			rank=rank-1;
			next:
			for(j=0;j<a.cols;j++){
				t[j]=a.data[r][j];
				a.data[r][j]=a.data[i][j];
				a.data[i][j];
			}
		}
    }
    return rank;
    return 0;
}

double trace_matrix(Matrix a)
{
    int i,tra;
    tra=0;
    if(a.rows==a.cols){
        for(i=0;i<a.rows;i++){
            tra=tra+a.data[i][i];
        }
        return tra;
    }else{
        printf("Error: The matrix must be a square matrix.");
    }
    return 0;
}

void print_matrix(Matrix a)
{
    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.cols; j++)
        {
            printf("%-8.2f", a.data[i][j]);
        }
        printf("\n");
    }
}
