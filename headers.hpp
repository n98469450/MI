#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <random>
#include <vector>
#include <cmath>
//#include <cstdint>
#define _USE_MATH_DEFINES

#define EPS 1e-300//я не хочу проверять, как при таком эпсилоне будут сходиться алгоритмы кластеризации...

using namespace std;


class Point {        

    private:

        double x;   
        double y;
        int number; 
        int colour; 

    public:

        Point();    

        double GetX();   
        double GetY();
        int GetColour();
        int GetNumber();

        void Set(double _x, double _y, int i, int c);
        void SetX(double _x);  
        void SetY(double _y);
        void SetNumber(int _n);
        void SetColour(int _c);

        ~Point();   

};


class Cloud {    

    private:

        int size;   
        Point center;    
        vector<Point> points; 
        vector<double> eigenvector1;
        vector<double> eigenvector2;

    public:

        Cloud(); 

        int GetSize(); 
        vector<Point> GetPoints();
        Point GetCenter();

        void AddPoint(Point AddPoint); //установка параметров облака (same)
        void AddCenter(Point c);
        void AddPoints(vector<Point> AddPoints);

        void Clear();    

};


class Cluster {   

    private:

        int size;
        Point center;
        vector<Point> points;
        vector<double> eigenvector1;
        vector<double> eigenvector2;

    public:

        Cluster(); 

        int GetSize();  
        vector<Point> GetPoints();
        Point GetCenter();
        vector<double> GetEig1();
        vector<double> GetEig2();

        void AddPoint(Point AddPoint); 
        void AddPoints(vector<Point> AddPoints);
        void AddCenter(Point c);
        void AddEig1(vector<double> eig_);
        void AddEig2(vector<double> eig_);

        void Clear(); //очищаем

};


class Field {  
    
    private:

	    int size;  
	    int condition; //информация о том, начата на поле кластеризация, или ещё нет
        int id;
        int IsSaved;
        int AmountOfCl;//количество проведённых на поле кластеризаций
        Cloud c;
    	vector<Cloud> clouds;
	    vector<Point> points;
        vector<vector<int>> incidence; 
        Point center;               
        vector<double> eigenvalues;
        vector<double> eigenvector1;
        vector<double> eigenvector2;

    public:
	
        Field();  

	    void AddCloud(Cloud AddCloud); 
	    void AddPoint(Point AddPoint);
	    void ReplaceCloudByNumber(int k, Cloud RepCloud);
	    void UpdatePoints();
        void SetTree(vector<vector<int>> in);
        void AddCenter(Point c);
        void AddEig1(vector<double> eig_);
        void AddEig2(vector<double> eig_);
        void SetID(int n);
        void SetIsSaved(int is);
        void SetCondition(int co);
        void SetAmountOfCl(int a);
        void SetEig(vector<double> eig);
    	
        
        int GetSize();        
        int GetCondition();
        int GetIsSaved();
        int GetID();
        int GetAmountOfCl();
        double distance (Point a, Point b);
        vector<vector<int>> GetTree();
        Point GetCenter();
	    Cloud GetCloudByNumber(int k);
	    vector<Point> GetPoints();
	    vector<Cloud> GetClouds();
        vector<double> GetEig1();
        vector<double> GetEig2();
        vector<double> GetEig();
    	
        void Clear();  

};


class Buffer {   
                 
    private:     
	    
        vector<Point> points;   
	    vector<double> center;  
        vector<double> eigenvector1; 
        vector<double> eigenvector2;
        vector<double> eigenvalues;

    public:

	    void AddCloud(Cloud AddCloud); 
        void AddPoints(vector<Point> AddPoints); 
        void FindFactors(vector<Point> points_); 
        void PrintFactors();                     
	    void Shift(double x, double y);  
	    void Rotate(double angle);       

        vector<double> GetCenter();   
        vector<double> GetEig1();
        vector<double> GetEig2();
        vector<double> GetEig();

        Cloud GetCloud();    
	
};


class Storage {   

    private:

        int size;  
        int IsSaved; //костыль, если возвращать как было -- закомментить   
	    vector<Cluster> clusters;
        int type;    
        int f;       
        double eps;
        vector<vector<int>> pseudo;  //НЕЙМИНГ //кстати, оно уже не нужно в коде, но выпиливать мне оч лень

    public:
	    
        Storage();  
	    
        void InitializeClusters(int k);  
	    void AddClusters(vector<Cluster> AddCluster); 
	    void AddCenterToCluster(int n, Point center);
	    void AddPointToCluster(int n, Point point);  
        void SetType(int t);
        void SetF(int f_);
        void SetEps(double eps_);
        void SetSize(int s);
        void SetSt(vector<vector<int>> st_);
        void UpdateClusters(Field field);
        void AddCluster(Cluster cluster_);
    	
        Point GetCenterFromCluster(int n);   
	    int GetSize();                       
    	vector<Cluster> GetClusters();        
        int GetType();
        int GetF();
        double GetEps();
        vector<vector<int>> GetSt();

        int GetIsSaved();//костыли
        void SetIsSaved(int IsSaved_);

	    void ClearAllClusters();   
    	void Clear();

};

//Далее -- алгоритмы кластеризации

class K_means {   

    private:

        Storage storage;

    public:

        void K_means_ (int k, Field field);
        Storage GetStorage();
};


class SpanningTree { 

    private:

        vector<vector<int>> incidence;

    public:

        void SpanningTree_(Field field);
        vector<vector<int>> GetIncidence();

};


class Wave {    

    private:

        Storage storage;
        Storage WaveFront;

    public: 
        
        void Wave_(double eps, Field field);
        Storage GetStorage();
        Storage GetWaveFront();

};


class EM {      

    private:

        Storage storage;     
        vector<vector<vector<double>>> cov;
        vector<Point> centroids;     

    public:

        void EM_(int k, Field field);
        Storage GetStorage();
        vector<vector<vector<double>>> GetCov(); 
        vector<Point> GetCentroids();            
        
};                                  


class DBSCAN {   //не работает нормально при многократном прогоне при одном запуске
                
    private:

        Storage storage;     
        Storage WaveClust;   

    public:

        void DBSCAN_(int k, double eps, Field field);
        Storage GetStorage();
        Storage GetWaveClust();

};


class Hierarchical {

    private:

        Storage storage;

        vector<vector<double>> distances;
        vector<Point> points;
        vector<int> ind;
        vector<vector<int>> sets;
        vector<vector<int>> structure;

    public:

        Storage GetStorage();
        
        void hierarchical(Field field, int parameter);

        double dist(Point x1, Point x2);
        
        void d_min(int step, int n);
        void d_max(int step, int n);
        void d_mid(int step, int n);
        void PrintTree();
        void Clusterise(int k);

};


class RBF{

    private:

        vector<double> ms;

    public:

        vector<double> GetMS();
        void SetMS(vector<double> ms_);
        
        double Regression(vector<Point> points, vector<double> meanings, double bandwidth, Point app);
        void Generate(int amount, double center, double disp);
        void Inter(vector<Point> points, vector<double> meanings, double bandwith, int sampling);
        double Regression1(vector<Point> points, double bandwidth, double x_);//функции под одномерную регрессию: её понятней, как рисовать
        void Inter1(vector<Point> points, double bandwidth, int sampling);

};


class Exec {    //хранит все сделанные за время кластеризации

    private:
            
        vector<Storage> storages; 

    public:
    
        K_means k_means;
        Wave wave;
        SpanningTree sp_tree;
        EM em;
        DBSCAN db;
        Hierarchical hc;
        RBF rbf;

        vector<Storage> GetStorages(); 
        void AddStorage(Storage st);
        void PrintEllipses(vector<vector<vector<double>>> cov, vector<Point> centers);
        void SetStorages(vector<Storage> st_);
        void REplaceStorage(int i, Storage st_);

        void SetSaved();

        void Clear();

};


class DataBase_Field {

    private:

        Field field;

    public:

        void SetField(Field field);
        Field GetField();

        void Save_Field(Field field);
        void UpdateAmount();
        int Read_Field(int number);
        int GetFieldNumber();           //при необходимости можно добавить стринговые аргументы, чтоб был выбор
};                                      //в названии файла-базы


class DataBase_Cluster {

    private:

        Cluster cluster;

    public:

        void SetCluster(Cluster cluster_);
        Cluster GetCluster();

        void Save_Cluster(Cluster cluster_, int FNumber, int ClNumber, int EtEN, int number, int TypeOf);//что за кластер сохраняется, номер поля, номер кластеризации
        void Read_Cluster(int FNumber, int ClNumber, int number);//номер поля, номер кластеризации, номер кластера
        int UpdateAmount();//возвращает количество сохранённых кластеров (надо для сквозной нумерации)

};


class DataBase_Storage {

    private:

        int ClN;//номер кластера при сквозной нумерации
        int EtEN;//номер кластеризации при сквозной нумерации
        int IsUpdated;
        int AmountOf;
        Storage storage;
        vector<Storage> storages;

    public:

        DataBase_Storage();
        
        DataBase_Cluster dbc;//возможно, стоит перенести в private (извне с ним всё равно не работают)

        void SetStorage(Storage storage_);
        Storage GetStorage();
        void SetStorages(vector<Storage> storages_);
        vector<Storage> GetStorages();
        void SetAmountOf(int AmountOf_);

        void Save_Storage(Storage storage_, int FNumber, int number);
        void Read_Storage(int FNumber, int number);
        void UpdateAmount();//возвращает число сохранённых FindCluster'ов (надо для сквозной нумерации)
        void Save_Storages(vector<Storage> storage_, int FNumber);
        void Read_Storages(int FNumber);
        
};


class Controller { 
                
    private:
	    
        int status;

    public:
	    
        Field field;
	    Storage storage;
    	Buffer buffer;
        Exec exec;
        DataBase_Field dbf;
        DataBase_Storage dbs;

    	Controller();
	                            
        vector<Point> Cloud_Gener(int n, double centerX, double centerY, double dispX, double dispY, int colour);
    	int Cloud_Find(int colour); 
        void SetField(Field field_);
	    void Cloud_Print(Cloud cloud);
    	void Field_Print();
	    void Print_to_File(int number);
        void Print_Tree();
	    
        ~Controller();

};


class Interface { 

    public:
	
        Controller controller;
	    
        int Start(string s);

};