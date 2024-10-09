#define _CRT_SECURE_NO_WARNINGS
#include "Header.h"
using namespace std;
typedef double type;

void input(string file, int N, type* mas)
{
    ifstream in;
    in.open(file + ".txt");
    for (int i = 0; i < N; i++)
        in >> mas[i];
    in.close();
}

void input_int(string file, int N, int* mas)
{
    ifstream in;
    in.open(file + ".txt");
    for (int i = 0; i < N; i++)
        in >> mas[i];
    if (mas[0])
        for (int i = 0; i <= N; i++)
            mas[i]--;
    in.close();
}

void output(string file, int N, type* mas)
{
    ofstream out;
    out.open(file + ".txt");
    for (int i = 0; i < N; i++)
        out << setprecision(17) << mas[i] << endl;
    out.close();
}



void LOS(int N, int maxit, type e, type* r, type* p, type* x, type* z, int* ia, int* ja, type* di, type* al, type* au, type* vec, type* Ar, type* res)
{
    type a, b, norm;
    int k = 0;
    type norm_v = sqrt(ScalarMultiplication(vec, vec, N));
    MatrixXVector(ia, ja, di, al, au, x, r, N);
    for (int i = 0; i < N; i++)
    {
        r[i] = vec[i] - r[i];
        z[i] = r[i];
    }
    MatrixXVector(ia, ja, di, al, au, z, p, N);
    for (k = 0; k < maxit && (sqrt(ScalarMultiplication(r, r, N)) / norm_v) > e; k++)
    {
        norm = ScalarMultiplication(p, p, N);
        a = ScalarMultiplication(p, r, N) / norm;
        CoeffXVector(z, a, res, N);
        VectorPlusVector(x, res, x, N);
        CoeffXVector(p, -a, res, N);
        VectorPlusVector(r, res, r, N);
        MatrixXVector(ia, ja, di, al, au, r, Ar, N);
        b = -ScalarMultiplication(p, Ar, N) / norm;
        CoeffXVector(z, b, res, N);
        VectorPlusVector(r, res, z, N);
        CoeffXVector(p, b, res, N);
        VectorPlusVector(Ar, res, p, N);
    }
    std::cout << "LOS" << "\n" << "Number of iterations: " << k << "\n" << "otnos nevyazka: " << (sqrt(ScalarMultiplication(r, r, N)) / norm_v) << "\n";
}


void LOS_D(int N, int maxit, type e, type* r, type* p, type* x, type* z, int* ia, int* ja, type* di, type* al, type* au, type* vec, type* LAUr, type* res, type* L)
{
    type a, b, norm;
    int m = 0;
    type norm_v = sqrt(ScalarMultiplication(vec, vec, N));
    for (int i = 0; i < N; i++)
        L[i] = 1 / sqrt(di[i]);
    MatrixXVector(ia, ja, di, al, au, x, res, N);
    for (int i = 0; i < N; i++)
        r[i] = L[i] * (vec[i] - res[i]);
    VectorXVector(L, r, z, N);
    MatrixXVector(ia, ja, di, al, au, z, res, N);
    VectorXVector(L, res, p, N);
    for (m = 0; m < maxit && (sqrt(ScalarMultiplication(r, r, N)) / norm_v) > e; m++)
    {
        norm = ScalarMultiplication(p, p, N);
        a = ScalarMultiplication(p, r, N) / norm;
        CoeffXVector(z, a, res, N);
        VectorPlusVector(x, res, x, N);
        CoeffXVector(p, -a, res, N);
        VectorPlusVector(r, res, r, N);
        VectorXVector(L, r, res, N);
        MatrixXVector(ia, ja, di, al, au, res, LAUr, N);
        VectorXVector(L, LAUr, LAUr, N);
        b = -ScalarMultiplication(p, LAUr, N) / norm;
        CoeffXVector(z, b, z, N);
        VectorXVector(L, r, res, N);
        VectorPlusVector(res, z, z, N);
        CoeffXVector(p, b, p, N);
        VectorPlusVector(p, LAUr, p, N);

    }
    std::cout << "LOSD" << "\n" << "Number of iterations: " << m << "\n" << "otnos nevyazka: " << (sqrt(ScalarMultiplication(r, r, N)) / norm_v) << "\n";
}

int Write_Txt_File_Of_Double(const char* fname, double* massiv, int n_of_records, int len_of_record)
{
    FILE* fp;
    int i, j;

    if ((fp = fopen(fname, "w")) == 0)
    {
        printf("Error: Cannot open file \"%s\" for writing.\n", fname);
        return 1;
    }

    printf("writing %s... ", fname);

    for (i = 0; i < n_of_records; i++)
    {
        for (j = 0; j < len_of_record; j++)
            fprintf(fp, "%25.13e\t", massiv[i * len_of_record + j]);
        fprintf(fp, "\n");
    }
    printf("done\n");

    fclose(fp);
    return 0;
}
//-------------------------------	
int Read_Long_From_Txt_File(const char* fname, int* number)
{
    FILE* fp;
    int temp;
    int retcode;

    if ((fp = fopen(fname, "r")) == 0)
    {
        char str[100];
        sprintf(str, "Error: Cannot open file \"%s\" for reading.\n", fname);
        return 1;
    }

    retcode = fscanf(fp, "%ld", &temp);
    if (retcode != 1)
    {
        char str[100];
        sprintf(str, "Error reading file \"%s\".\n", fname);
        fclose(fp);
    }
    *number = temp;

    fclose(fp);
    return 0;
}
//-------------------------------
int Read_Bin_File_Of_Double(const char* fname, double* massiv, int n_of_records, int len_of_record)
{
    int temp;
    FILE* fp;

    // открываем файл на чтение
    if ((fp = fopen(fname, "r+b")) == 0)
    {
        char str[100];
        sprintf(str, "Error: Cannot open file \"%s\" for reading.\n", fname);
        return 1;
    }

    temp = fread(massiv, sizeof(double) * len_of_record, n_of_records, fp);
    if (temp != n_of_records)
    {
        char str[100];
        sprintf(str, "Error reading file \"%s\". %ld of %ld records was read.\n", fname, temp, n_of_records);
        fclose(fp);
        return 1;
    }

    fclose(fp);

    return 0;
}
//-------------------------------
int Read_Bin_File_Of_Long(const char* fname, int* massiv, int n_of_records, int len_of_record)
{
    int temp;
    FILE* fp;

    // открываем файл на чтение
    if ((fp = fopen(fname, "r+b")) == 0)
    {
        char str[100];
        sprintf(str, "Cannot open file %s.\n", fname);
        return 1;
    }

    // чтение
    temp = fread(massiv, sizeof(int) * len_of_record, n_of_records, fp);
    if (temp != n_of_records)
    {
        char str[100];
        sprintf(str, "Error reading file \"%s\". %ld of %ld records was read.\n", fname, temp, n_of_records);
        fclose(fp);
        return 1;
    }


    fclose(fp);
    return 0;
}

//------------------------------------------------------------
void FromRSFToCSR_Real_1_Sym(int nb, int* ig, int* sz_ia, int* sz_ja)
{
    *sz_ia = nb + 1;
    *sz_ja = ig[nb] + nb;
}
//------------------------------------------------------------------------

void main()
{
    /*ifstream in;
    in.open("kuslau.txt");
    int N, maxit;
    type e;
    in >> N >> maxit >> e;
    in.close();

    type* x = new type[N];
    for (int i = 0; i < N; i++)
        x[i] = 0;

    int* ig = new int[N + 1];
    input_int("ig.txt", N + 1, ig);

    int size = ig[N];
    int* jg = new int[size];
    type* ggl = new type[size], * ggu = new type[size];
    type* di = new type[N], * pr = new type[N], * res = new type[N];
    type* L = new type[N], * Ar = new type[N];
    type* r = new type[N], * z = new type[N], * p = new type[N];
    input("di.txt", N, di);
    input_int("jg.txt", size, jg);
    input("gg.txt", size, ggl);
    input("gg.txt", size, ggu);
    input("pr.txt", N, pr);*/


    // ######################### FOR BINARY FILES #################################


    int N = 0;
    int maxit = 0;
    type e = 0;
    // nb
    ifstream in;
    in.open("kuslau.txt");
    in >> N >> maxit >> e;
    in.close();
    int* ig = new int[N + 1];
    int ig_n_1 = 0;
    Read_Bin_File_Of_Long("ig", ig, N + 1, 1);
    for (int i = 0; i < N + 1; i++)
    {
        ig[i]--;
    }
    ig_n_1 = ig[N];
    int size = ig[N];
    int* jg = new int[size];
    type* ggl = new type[size], * ggu = new type[size];
    type* di = new type[N], * pr = new type[N], * res = new type[N];
    type* L = new type[N], * Ar = new type[N];
    type* r = new type[N], * z = new type[N], * p = new type[N];
    type* x = new type[N];
    for (int i = 0; i < N; i++)
        x[i] = 0;


    ////////////////////////////////////////
    int sz_ia = 0;
    int sz_ja = 0;


    // jg
    jg = new int[ig_n_1];
    Read_Bin_File_Of_Long("jg", jg, ig_n_1, 1);
    for (int i = 0; i < ig_n_1; i++)
    {
        jg[i]--;
    }

    // di
    di = new double[N];
    Read_Bin_File_Of_Double("di", di, N, 1);

    // ggl
    ggl = new double[ig_n_1];
    Read_Bin_File_Of_Double("gg", ggl, ig_n_1, 1);

    ggu = new double[ig_n_1];
    Read_Bin_File_Of_Double("gg", ggu, ig_n_1, 1);



    // pr
    pr = new double[N];
    Read_Bin_File_Of_Double("pr", pr, N, 1);

    std::cout << "1-LOS\n" << "2-LOS_D\n" << "3-PARDISO\n";
    int a = 0;
    std::cin >> a;
    switch (a)
    {
    case 1:
    {
        clock_t start = clock();
        LOS(N, maxit, e, r, p, x, z, ig, jg, di, ggl, ggu, pr, Ar, res);
        output("out", N, x);
        clock_t end = clock();
        double seconds = (double)(end - start) / CLOCKS_PER_SEC;
        printf("The time: %f seconds\n", seconds);
        break;
    }
    case 2:
    {
        clock_t start = clock();
        LOS_D(N, maxit, e, r, p, x, z, ig, jg, di, ggl, ggu, pr, Ar, res, L);
        output("out", N, x);
        clock_t end = clock();
        double seconds = (double)(end - start) / CLOCKS_PER_SEC;
        printf("The time: %f seconds\n", seconds);
        break;
    }
    default:
        break;
    }
}
