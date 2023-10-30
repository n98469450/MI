#include "headers.hpp"


double d(Point a, Point b){

    double aX, bX;
    double aY, bY;
    double distance = 0.0;

    aX = a.GetX(); bX = b.GetX();
    aY = a.GetY(); bY = b.GetY();

    distance = sqrt (pow(aX - bX, 2) + pow(aY - bY, 2));

    return distance;
}


//Point namespace

Point::Point(){

    x = 0;
    y = 0;
    number = 0;
    colour = 0;

}

double Point::GetX(){ return x; }
double Point::GetY(){ return y; }
int Point::GetColour(){ return colour; }
int Point::GetNumber(){ return number; }

void Point::Set(double _x, double _y, int i, int c){ x = _x; y = _y; number = i; colour = c; }
void Point::SetX(double _x){ x = _x; }
void Point::SetY(double _y){ y = _y; }
void Point::SetNumber(int _n){ number = _n; }
void Point::SetColour(int _c){ colour = _c; }

Point::~Point(){}


//Cloud namespace

Cloud::Cloud(){ size = 0; }

void Cloud::AddPoint(Point AddPoint){
        
    points.push_back(AddPoint);
    size = points.size();

}

void Cloud::AddPoints(vector<Point> AddPoints){
    
    points = AddPoints;
    size = points.size();  

}

void Cloud::AddCenter(Point c){ center = c; }

int Cloud::GetSize(){ return size; }
vector<Point> Cloud::GetPoints(){ return points; }
Point Cloud::GetCenter(){ return center; }

void Cloud::Clear(){ points.clear(); }


//Cluster namespace

Cluster::Cluster(){

    size = 0;
    center.Set(0, 0, 0, 0);

}

void Cluster::AddCenter(Point AddCenter){ center = AddCenter; }

void Cluster::AddPoint(Point AddPoint){

    points.push_back(AddPoint);
    size = points.size();

}

void Cluster::AddPoints(vector<Point> AddPoints){

    points.insert(points.end(), AddPoints.begin(), AddPoints.end());
    size = points.size();

}

void Cluster::AddEig1(vector<double> eig_){ eigenvector1 = eig_; }
//void Cluster::AddEig2(vector<double> eig_){ eigenvector1 = eig_; }//а, ясно
                                                        //это стирать нельзя, останется памятником моему [ДАННЫЕ УДАЛЕНЫ]
void Cluster::AddEig2(vector<double> eig_){ eigenvector2 = eig_; }
int Cluster::GetSize(){ return size; }

vector<double> Cluster::GetEig1(){ return eigenvector1; }
vector<double> Cluster::GetEig2(){ return eigenvector2; }
Point Cluster::GetCenter(){ return center; }
vector<Point> Cluster::GetPoints(){ return points; }

void Cluster::Clear(){

    size = 0;
    points.clear();

}


//Field namespace

Field::Field(){

    size = 0;
    condition = 0;
    id = -1;
    IsSaved = 0;
    AmountOfCl = 0;

}

void Field::SetID(int n){ id = n; }
void Field::SetIsSaved(int is){ IsSaved = is; }
void Field::SetCondition(int co){ condition = co; }
void Field::SetAmountOfCl(int a){ AmountOfCl = a; }

void Field::AddCloud(Cloud AddCloud){

    vector<Point> AddPoints = AddCloud.GetPoints();
    clouds.push_back(AddCloud);
    points.insert(points.end(), AddPoints.begin(), AddPoints.end());
    size = points.size();

}

void Field::AddCenter(Point c){ center = c; }

void Field::AddPoint(Point AddPoint){ points.push_back(AddPoint); }

void Field::AddEig1(vector<double> eig_){ eigenvector1 = eig_; }
void Field::AddEig2(vector<double> eig_){ eigenvector2 = eig_; }
void Field::SetEig(vector<double> eig){ eigenvalues = eig; }

void Field::ReplaceCloudByNumber(int k, Cloud RepCloud){ clouds[k].AddPoints(RepCloud.GetPoints()); }

void Field::UpdatePoints(){

    vector<Point> UpdPoints;
    int n = clouds.size();
    int i, j, l = 0;

    points.resize(0);

    for (i = 0; i < n; i++){

        UpdPoints = clouds[i].GetPoints();

        for (j = 0; j < int(UpdPoints.size()); j++){

            points.push_back(UpdPoints[j]);
            l++;

        }
    }
}

int Field::GetSize(){ return points.size(); }

int Field::GetIsSaved(){ return IsSaved; }

int Field::GetAmountOfCl(){ return AmountOfCl; }

void Field::SetTree(vector<vector<int>> in){ incidence = in; }

int Field::GetCondition(){ return condition; }

Point Field::GetCenter(){ return center; }
vector<double> Field::GetEig1(){ return eigenvector1; }
vector<double> Field::GetEig2(){ return eigenvector2; }
vector<double> Field::GetEig(){ return eigenvalues; }

vector<vector<int>> Field::GetTree(){ return incidence; }

Cloud Field::GetCloudByNumber(int k){

    int n = clouds.size();
    Cloud cl;
    if (n - 1 < k) return cl;
    return clouds[k];

}

vector<Point> Field::GetPoints(){ return points; }

vector<Cloud> Field::GetClouds(){ return clouds; }

double Field::distance (Point a, Point b){
    
    double ax = a.GetX(), ay = a.GetY();
    double bx = b.GetX(), by = b.GetY();

    return sqrt(pow(ax - bx, 2) + pow (ay - by, 2));
}

int Field::GetID(){ return id; }

void Field::Clear(){

    clouds.clear();
    points.clear();
    
}


//Buffer namespace

void  Buffer::AddCloud(Cloud AddCloud){

    points = AddCloud.GetPoints();
    center[0] = AddCloud.GetCenter().GetX();
    center[1] = AddCloud.GetCenter().GetY();

}

void  Buffer::Shift(double x, double y){

    double p_x, p_y;
    int i;

    for (i = 0; i < int(points.size()); i++){

        p_x = points[i].GetX();
        p_y = points[i].GetY();

        p_x += x;
        p_y += y;
        points[i].Set(p_x, p_y, i, i);

    }
}

void  Buffer::Rotate(double angle){

    double p_x, p_y, new_x, new_y;
    int i;
    int colour;

    for (i = 0; i < int(points.size()); i++){

        p_x = points[i].GetX() - center[0];
        p_y = points[i].GetY() - center[1];

        new_x = cos(angle) * p_x - sin(angle) * p_y + center[0];
        new_y = sin(angle) * p_x + cos(angle) * p_y + center[1];
        colour = points[i].GetColour();
        points[i].Set(new_x, new_y, i, colour);

    }
}

Cloud Buffer::GetCloud(){

    Cloud cl;
    cl.AddPoints(points);

    return cl;
}

vector<double> Buffer::GetCenter(){ return center; }

void Buffer::AddPoints(vector<Point> AddPoints){ points = AddPoints; }

void Buffer::FindFactors(vector<Point> points_){
     
    vector<double> x;
    vector<double> y;
    double cx = 0.0, cy = 0.0;
    double a11 = 0.0, a12 = 0.0, a22 = 0.0;
    double l1, l2;
    int i, n;
    n = points_.size();

    x.resize(n);
    y.resize(n);
    eigenvalues.resize(0);

    for (i = 0; i < n; i++){
        cx += points_[i].GetX();
        cy += points_[i].GetY();
    }
    
    cx /= n;
    cy /= n;

    for (i = 0; i < n; i++){

        x[i] = points_[i].GetX() + cx;
        y[i] = points_[i].GetY() + cy;
    
    }

    for (i = 0; i < n; i++){

        a11 += x[i] * x[i];
        a12 += x[i] * y[i];
        a22 += y[i] * y[i];
        
    }

    l1 = (a11 + a22 + sqrt((a11 - a22) * (a11 - a22) + 4 * a12 * a12)) / 2;
    l2 = (a11 + a22 - sqrt((a11 - a22) * (a11 - a22) + 4 * a12 * a12)) / 2;

    eigenvalues.push_back(l1);
    eigenvalues.push_back(l2);

    center.resize(2);
    eigenvector1.resize(2);
    eigenvector2.resize(2);

    center[0] = cx;
    center[1] = cy;

    eigenvector1[0] = 1;//сюда ещё нормировку надо добавить на самом деле
    eigenvector1[1] = -(a11 - l1) / a12;
    eigenvector2[0] = 1;
    eigenvector2[1] = -(a11 - l2) / a12;

}

void Buffer::PrintFactors(){

    ofstream file;

    file.open("factors.txt");

    file << center[0] << " " << center[1] << endl;
    file << eigenvector1[0] << " " << eigenvector1[1] << endl;
    file << eigenvector2[0] << " " << eigenvector2[1] << endl;
    file << endl;

    file.close();
}

vector<double> Buffer::GetEig1(){ return eigenvector1; }
vector<double> Buffer::GetEig2(){ return eigenvector2; }
vector<double> Buffer::GetEig(){ return eigenvalues; }


//Storage namespace

Storage::Storage(){ size = 0; IsSaved = 0; }

void Storage::InitializeClusters(int k){

    clusters.resize(k);
    size = k;

}

void Storage::AddClusters(vector<Cluster> AddClusters){ //интересный факт: эта функция не используется НИГДЕ

    clusters = AddClusters;
    size = AddClusters.size();

}

void Storage::AddCenterToCluster(int n, Point center){ clusters[n].AddCenter(center); }
void Storage::AddPointToCluster(int n, Point point){ clusters[n].AddPoint(point); }
Point Storage::GetCenterFromCluster(int n){ return clusters[n].GetCenter(); }
int Storage::GetSize(){ return size; }
vector<Cluster> Storage::GetClusters(){ return clusters; }
vector<vector<int>> Storage::GetSt(){ return pseudo; }     //НЕЙМИНГ
void Storage::ClearAllClusters(){ for (int i = 0; i < size; i++){ clusters[i].Clear(); } IsSaved = 0; }
int Storage::GetType(){ return type; }
int Storage::GetF(){ return f; }
double Storage::GetEps(){ return eps; }
void Storage::SetType(int t){ type = t; }
void Storage::SetF(int f_){ f = f_; }
void Storage::SetEps(double eps_){ eps = eps_; }
void Storage::SetSize(int s){ size = s; }
void Storage::SetSt(vector<vector<int>> st_){ pseudo = st_; }
void Storage::AddCluster(Cluster cluster_){ clusters.push_back(cluster_); }

int Storage::GetIsSaved(){ return IsSaved; };
void Storage::SetIsSaved(int IsSaved_){ IsSaved = IsSaved_; }

void Storage::UpdateClusters(Field field){

    int i, j, n;
    int size = pseudo.size();

    clusters.clear(); //для надёжности
    InitializeClusters(size);

    for (i = 0; i < size; i++){

        n = pseudo[i].size();

        for (j = 0; j < n; j++){ AddPointToCluster(i, field.GetPoints()[pseudo[i][j]]); }

    }
}

void Storage::Clear(){ clusters.clear(); IsSaved = 0; }


//Controller namespace

Controller::Controller(){ status = 1; }

vector<Point> Controller::Cloud_Gener(int n, double centerX, double centerY, double dispX, double dispY, int colour){
        
    vector<Point> points;
    Point a;
    default_random_engine gener;
    normal_distribution<double> distributionx(centerX, dispX);
    normal_distribution<double> distributiony(centerY, dispY);

    for (int i = 0; i < n; i++){

        a.Set(distributionx(gener), distributiony(gener), i, colour);
        points.push_back(a);

    }

    return points;
}

void Controller::SetField(Field field_){ field = field_; }

void Controller::Cloud_Print(Cloud cloud){
        
    int n = cloud.GetSize();
    vector<Point> arr = cloud.GetPoints();

    for (int i = 0; i < n; i++){
        cout << arr[i].GetX() << " " << arr[i].GetY() << endl;
    }

    cout << endl;

    return;
}

void Controller::Field_Print(){

    vector<Point> point = field.GetPoints();
    int n = point.size();
    ofstream file;

    file.open("field.txt");

    for (int i = 0; i < n; i++){
        file << point[i].GetX() << " " << point[i].GetY() << endl;
    }

    file.close();

    return;
}

void Controller::Print_to_File(int number){

    ofstream file1, file2;
    vector<Cluster> clusters = exec.GetStorages()[number].GetClusters();
    vector<Point> points;
    stringstream num;
    string pref;
    string pref_c = "_";
    string index = "\0";
    string post = ".txt";
    string s;
    int k = clusters.size();
    int i, j, n;

    num.str().resize(0);
    num.str().clear();
    num.str("");
    num << (number + 1);
    pref = num.str();

    file2.open("factors.txt");

    for (i = 0; i < k; i++){

        num.str().resize(0);
        num.str().clear();
        num.str("");
        num << (i + 1);
        index = num.str();
        s = pref + pref_c + index + post;
        file1.open(s);

        points = clusters[i].GetPoints();
        n = clusters[i].GetSize();

        buffer.FindFactors(points);

        file2 << "set arrow from " << buffer.GetCenter()[0] << "," << buffer.GetCenter()[1] << " to ";// << endl;
        file2 << buffer.GetEig1()[0] + buffer.GetCenter()[0] << "," << buffer.GetEig1()[1] + buffer.GetCenter()[1] << " lw 2" << endl;
        file2 << "plot (0,0) lc 0" << endl;
        file2 << "set arrow from " << buffer.GetCenter()[0] << "," << buffer.GetCenter()[1] << " to ";// << endl;
        file2 << buffer.GetEig2()[0] + buffer.GetCenter()[0] << "," << buffer.GetEig2()[1] + buffer.GetCenter()[1] << " lw 2" << endl;
        file2 << "plot (0,0) lc 0" << endl;
        file2 << endl;

        for (j = 0; j < n; j++){
            file1 << points[j].GetX() << " " << points[j].GetY() << endl;
        }

    file1.close();

    }
    file2.close();
}

void Controller::Print_Tree(){

    ofstream file;
    ofstream l;
    vector<Point> points = field.GetPoints();
    vector<vector<int>> tree = field.GetTree();
    int n = points.size();
    int i, j;

    file.open("tree.txt");
    l.open("data.txt");


    for (i = 0; i < n; i++){
        for (j = i + 1; j < n; j++){
            if (tree[i][j]){
                file << points[i].GetX() << " " << points[i].GetY() << endl;
                file << points[j].GetX() << " " << points[j].GetY() << endl;
                file << endl;
                l << sqrt((points[i].GetX() - points[j].GetX())*(points[i].GetX() - points[j].GetX()) + (points[i].GetX() - points[j].GetX())*(points[i].GetX() - points[j].GetX()));
                l << endl;
            }
        }
    }

    file.close();

}

Controller::~Controller(){ }//загадочные обстоятельства


//Interface namespace

int	Interface::Start(string s){

    string s0 = "help";
    string s1 = "create_cloud";
    string s2 = "print_cloud";
    string s3 = "print_field";
    string s4 = "k-means";
    string s5 = "print_clusters_to_files";
    string s6 = "wave";
    string s7 = "dbscan";
    string s8 = "rotate";
    string s9 = "stop";
    string s10 = "spanning_tree";
    string s11 = "print_tree";
    string s12 = "em";
    string s13 = "print_ellipses";
    string s14 = "print_factors";
    string s15 = "info";
    string s16 = "save_field";
    string s17 = "clear_field";
    string s18 = "read_field";
    string s19 = "save_storages";   //надо чёт придумать с этими командами, а то чёт путаная немного терминология выходит
    string s20 = "clear_storages";
    string s21 = "read_storages";
    string s22 = "hierarchical";
    string s23 = "read_storage";
    string s24 = "regression";
    string s25 = "interpolation";
    string s26 = "regression_1d";
    string s27 = "interpolation_1d";

    if (s == s0) return 0;
    if (s == s1) return 1;
    if (s == s2) return 2;
    if (s == s3) return 3;
    if (s == s4) return 4;
    if (s == s5) return 5;
    if (s == s6) return 6;
    if (s == s7) return 7;
    if (s == s8) return 8;
    if (s == s9) return 9;
    if (s == s10) return 10;
    if (s == s11) return 11;
    if (s == s12) return 12;
    if (s == s13) return 13;
    if (s == s14) return 14;
    if (s == s15) return 15;
    if (s == s16) return 16;
    if (s == s17) return 17;
    if (s == s18) return 18;
    if (s == s19) return 19;
    if (s == s20) return 20;
    if (s == s21) return 21;
    if (s == s22) return 22;
    if (s == s23) return 23;
    if (s == s24) return 24;
    if (s == s25) return 25;
    if (s == s26) return 26;
    if (s == s27) return 27;
    

    return -1;
}


//K_means namespace

Storage K_means::GetStorage(){ return storage; }

void K_means::K_means_ (int k, Field field){

    vector<vector<int>> ind;
    vector<double> dist;
    vector<Point> points = field.GetPoints();
    vector<Point> centroids1, centroids2;
    double MinDist;
    double sumX, sumY;
    int amount; 
    int i, j, c, flag = 1;
    int n = points.size();

    //cout << n << endl;

    ind.resize(k);
    for (i = 0; i < k; i++){
        ind[i].resize(n);
    }
    centroids1.resize(k);
    centroids2.resize(k);
    dist.resize(k);

    for (i = 0; i < k; i++){
        centroids1[i] = points[i];
    }

    while (flag){
        for (i = 0; i < n; i++){
            for (j = 0; j < k; j++){
                dist[j] = field.distance(points[i], centroids1[j]);
            }
            MinDist = dist[0]; c = 0;
            for (j = 1; j < k; j++){
                if (MinDist > dist[j]){ MinDist = dist[j]; c = j; }
            }

            for (j = 0; j < k; j++){
                if (fabs(MinDist - dist[j]) < EPS){ ind[j][i] = 1; }
                else { ind [j][i] = 0; }
            }
        }
        
        for (i = 0; i < k; i++){
            amount = 0;
            sumX = 0;
            sumY = 0;
            for (j = 0; j < n; j++){
                sumX += ind[i][j]*points[j].GetX();
                sumY += ind[i][j]*points[j].GetY();
                amount += ind[i][j];
            }
            centroids2[i].SetX(sumX/amount);
            centroids2[i].SetY(sumY/amount);
        }
        flag = 0;
        for (i = 0; i < k; i++){
            if (field.distance(centroids1[i], centroids2[i]) > EPS){ flag = 1; }
            centroids1[i] = centroids2[i];
        }
    }

    //int correct = 0;
    //for (i = 0; i < k; i++){ for (j = 0; j < n; j++){ correct += ind[i][j]; } cout << correct << "\t"; correct = 0; }
    //cout << endl << endl;

    storage.Clear();
    
    storage.InitializeClusters(k);
    for (i = 0; i < k; i++){
        for (j = 0; j < n; j++){
            if (ind[i][j] == 1){
                storage.AddPointToCluster(i, points[j]);
            }
        }
    }
    storage.SetType(1);
}


//Wave namespace

Storage Wave::GetStorage(){ return storage; }

Storage Wave::GetWaveFront(){ return WaveFront; }

void Wave::Wave_(double eps, Field field){

    vector<Point> points = field.GetPoints();
    vector<vector<double>> dist;
    vector<int> wave_front; 
    vector<int> ind;        
    int n = points.size();  
    int i, j;               
    int iteration = 1;      
    int flag = 1;           
    int AmountOfClusters = 0;

    dist.resize(n);
    for (i = 0; i < n; i++){ dist[i].resize(n); }
    wave_front.resize(n);
    ind.push_back(0);

    
    for (i = 0; i < n; i++){  
        for (j = 0; j < n; j++){ dist[i][j] = field.distance(points[i], points[j]); }
    }

    for (i = 0; i < n; i++){  
        for (j = 0; j < n; j++){ 
            if (dist[i][j] <= eps){ dist [i][j] = 1; }
            else{ dist[i][j] = 0; }
        }
    }

    while(1){

        for (i = 1; i < n; i++){
            if (wave_front[i] == 0){

                wave_front[i] = iteration;
                flag = 1;
                break;

            }
        }

        if (flag == 0){ break; }

        while (flag){  

            flag = 0;

            for (i = 0; i < n; i++){
                if (wave_front[i] == iteration){
                    for (j = 0; j < n; j++){
                        if (dist[i][j] && (wave_front[j] == 0)){
                            wave_front[j] = iteration + 1;
                            flag = 1;
                        }
                    }
                }
            }
            iteration++;
        }
        AmountOfClusters++;
        ind.push_back(iteration);
    }

    WaveFront.InitializeClusters(iteration);

    for (i = 1; i < iteration; i++){
        for (j = 0; j < n; j++){
            if (wave_front[j] == i){
                WaveFront.AddPointToCluster((i - 1), points[j]);
            }
        }
    }

    WaveFront.SetType(2);
    WaveFront.SetEps(eps);

    storage.InitializeClusters(AmountOfClusters);

    for (i = 0; i < n; i++){
        for (j = 0; j < AmountOfClusters; j++){
            if ((wave_front[i] >= ind[j]) && (wave_front[i] < ind[j+1])){
                storage.AddPointToCluster(j, points[i]);
            }
        }
    }

    storage.SetType(3);
    storage.SetEps(eps);
}


//Spanning tree namespace

vector<vector<int>> SpanningTree::GetIncidence(){ return incidence; }

void SpanningTree::SpanningTree_(Field field){

    vector<Point> points = field.GetPoints();
    vector<vector<double>> dist;
    vector<int> ind;
    vector<int> edge;
    double mindist;
    int n = points.size();
    int i, j;
    int flag = 1, first = 1;
    int iter = 0;

    dist.resize(n);
    for (i = 0; i < n; i++){ dist[i].resize(n); }
    incidence.resize(n);
    for (i = 0; i < n; i++){ incidence[i].resize(n); }
    edge.resize(2);
    ind.resize(n);

    for (i = 0; i < n; i++){ 
        for (j = 0; j < n; j++){ dist[i][j] = field.distance(points[i], points[j]); }
    }

    mindist = dist[0][1];
    edge[0] = 0;
    edge[1] = 1;

    for (i = 0; i < n; i++){
        for (j = i + 1; j < n; j++){
            if (dist[i][j] < mindist){
                mindist = dist[i][j];
                edge[0] = i;
                edge[1] = j;
            }
        }
    }

    ind[edge[0]] = 1;
    ind[edge[1]] = 1;
    incidence[edge[0]][edge[1]] = 1;
    incidence[edge[1]][edge[0]] = 1;

    while (flag){

        flag = 0;

        for (i = 0; i < n; i++){
            for (j = i + 1; j < n; j++){
                if (((ind[i] + ind[j]) > 0) && ((ind[i]*ind[j]) == 0)){
                    if (first){
                        mindist = dist[i][j];
                        edge[0] = i;
                        edge[1] = j;
                        first = 0;
                    }
                    else if (mindist > dist[i][j]){
                        mindist = dist[i][j];
                        edge[0] = i;
                        edge[1] = j;
                    }
                }
            }
        }

        ind[edge[0]] = 1;
        ind[edge[1]] = 1;
        incidence[edge[0]][edge[1]] = 1;
        incidence[edge[1]][edge[0]] = 1;

        for (i = 0; i < n; i++){
            if (ind[i] == 0){
                flag = 1;
                first = 1;
                break;
            }
        }
        iter++;
    }
}


//EM namespace

Storage EM::GetStorage(){ return storage; }
vector<vector<vector<double>>> EM::GetCov(){ return cov; }
vector<Point> EM::GetCentroids(){ return centroids; }

void EM::EM_(int k, Field field){

    vector<Point> points = field.GetPoints();
    vector<vector<double>> prob;
    vector<vector<double>> prev_prob;
    vector<double> P;
    vector<double> v;
    int n = points.size();
    int i, j;
    int flag = 1;
    int index;
    double MaxProb;
    double sigma12, sigma22;
    double sumP;
    double exponent;
    double detCov;
    double normalization;
    double wSumX, wSumY;
    double check;

    int iter = 0;

    prev_prob.resize(k);

    for (i = 0; i < k; i++){

        for (j = 0; j < n; j++){
            prev_prob[i].push_back(1.0/k);
        }

    }

    prob.resize(k);

    for (i = 0; i < k; i++){ prob[i].resize(n); }
    for (i = 0; i < k; i++){ P.push_back(1.0/k); }

    v.resize(2);
    
    for(i = 0; i < n; i++){

        sigma12 += (points[i].GetX()*points[i].GetX())/n;
        sigma22 += (points[i].GetY()*points[i].GetY())/n;

    }
    cov.resize(k);
    for (i = 0; i < k; i++){

        cov[i].resize(2);   

        cov[i][0].push_back(sigma12);
        cov[i][0].push_back(0);
        cov[i][1].push_back(0);
        cov[i][1].push_back(sigma22);

    }

    for (i = 0; i < k; i++){ centroids.push_back(points[i]); }

    while (flag){

        flag = 0;

        for (i = 0; i < k; i++){ 

            for (j = 0; j < n; j++){

                detCov = cov[i][0][0]*cov[i][1][1] - cov[i][0][1]*cov[i][1][0];
                v[0] = points[j].GetX() - centroids[i].GetX();
                v[1] = points[j].GetY() - centroids[i].GetY();
                exponent = v[0]*v[0]*cov[i][1][1] - v[0]*v[1]*(cov[i][0][1] + cov[i][1][0]) + v[1]*v[1]*cov[i][0][0];
                exponent /= (-2)*detCov;
                prob[i][j] = exp(exponent)/sqrt(2*M_PI*detCov);

            }
        }

        for (i = 0; i < n; i++){ //нормируем вероятности, чтоб сумма по строке была 1

            normalization = 0;

            for (j = 0; j < k; j++){ normalization += prob[j][i]*P[j]; }

            for (j = 0; j < k; j++){

                prob[j][i] /= normalization;
                prob[j][i] *= P[j];

            }
        }

        for (i = 0; i < k; i++){ 

            sumP = 0;

            for (j = 0; j < n; j++){ sumP += prob[i][j]; }

            P[i] = sumP/n;

        }

        
        for (i = 0; i < k; i++){

            wSumX = 0;
            wSumY = 0;

            for (j = 0; j < n; j++){

                wSumX += (prob[i][j]/P[i])*points[j].GetX();
                wSumY += (prob[i][j]/P[i])*points[j].GetY();

            }

            centroids[i].SetX(wSumX/n);
            centroids[i].SetY(wSumY/n);

        }
        //матрицы ковариации
        for (i = 0; i < k; i++){

            cov[i][0][0] = 0;
            cov[i][0][1] = 0;
            cov[i][1][0] = 0;
            cov[i][1][1] = 0;
            
            for (j = 0; j < n; j++){

                cov[i][0][0] += (prob[i][j]/P[i])*(points[j].GetX() - centroids[i].GetX())*(points[j].GetX() - centroids[i].GetX());
                cov[i][0][1] += (prob[i][j]/P[i])*(points[j].GetX() - centroids[i].GetX())*(points[j].GetY() - centroids[i].GetY());
                cov[i][1][0] += (prob[i][j]/P[i])*(points[j].GetY() - centroids[i].GetY())*(points[j].GetX() - centroids[i].GetX());
                cov[i][1][1] += (prob[i][j]/P[i])*(points[j].GetY() - centroids[i].GetY())*(points[j].GetY() - centroids[i].GetY());

            }

            cov[i][0][0] /= n;
            cov[i][0][1] /= n;
            cov[i][1][0] /= n;
            cov[i][1][1] /= n;

        }

        //условие схождения алгоритма
        for (i = 0; i < k; i++){

            check = 0;

            for (j = 0; j < n; j++){

                check += fabs(prob[i][j] - prev_prob[i][j]);
                prev_prob[i][j] = prob[i][j];

            }

            check /= n;

            if (check >= EPS){

                flag = 1;
                break;

            }
        }
    }

    storage.InitializeClusters(k);

    for (i = 0; i < n; i++){

        MaxProb = 0.0;
        index = 0;

        for (j = 0; j < k; j++){

            if (prob[j][i] > MaxProb){

                MaxProb = prob[j][i];
                index = j;

            }
        }
        storage.AddPointToCluster(index, points[i]);
    }
    storage.SetType(6);

}


//DBScan namespace

Storage DBSCAN::GetStorage(){ return storage; }
Storage DBSCAN::GetWaveClust(){ return WaveClust; }

void DBSCAN::DBSCAN_(int k, double eps, Field field){

    vector<Point> points = field.GetPoints();
    vector<vector<double>> dist;
    vector<int> ind;
    vector<int> index;
    vector<int> wave_front;
    int i, j;
    int iteration = 1;
    int flag = 1;
    int AmountOfClusters = 0;
    double f;
    int n = points.size();

    dist.resize(n);
    for (i = 0; i < n; i++){ dist[i].resize(n); }
    ind.resize(0);
    for (i = 0; i < n; i++){ ind.push_back(0); }
    index.clear();
    index.push_back(0);
    wave_front.resize(n);

    for (i = 0; i < n; i++){  
        for (j = 0; j < n; j++){ dist[i][j] = field.distance(points[i], points[j]); }
    }

    for (i = 0; i < n; i++){  
        for (j = 0; j < n; j++){ 
            if (dist[i][j] <= eps){ dist [i][j] = 1; }
            else{ dist[i][j] = 0; }
        }
    }

    for (i = 0; i < n; i++){  
        f = 0;
        for (j = 0; j < n; j++){
            f += dist[i][j];
        }
        if (f >= k){
            ind[i] = 1;
        }
    }

    for (i = 0; i < n; i++){  
        if (ind[i] == 0){
            for (j = 0; j < n; j++){
                if ((dist[i][j] == 1) && (ind[j] == 1)){
                    ind[i] = 2;
                    break;
                }
            }
        }
    }

    storage.InitializeClusters(3);

    for (i = 0; i < n; i++){
        if (ind[i] == 0){  
            storage.AddPointToCluster(2, points[i]);
        }
        else if (ind[i] == 1){
            storage.AddPointToCluster(0, points[i]);
        }
        else if (ind[i] == 2){
            storage.AddPointToCluster(1, points[i]);
        }
    }

    storage.SetType(4);
    storage.SetEps(eps);
    storage.SetF(k);

    WaveClust.SetType(5);
    WaveClust.SetEps(eps);
    WaveClust.SetF(k);

    while(1){

        for (i = 1; i < n; i++){
            if ((wave_front[i] == 0) && (ind[i] != 0)){

                wave_front[i] = iteration;
                flag = 1;
                break;

            }
        }

        if (flag == 0){ break; }

        while (flag){  

            flag = 0;

            for (i = 0; i < n; i++){
                if (wave_front[i] == iteration){
                    for (j = 0; j < n; j++){
                        if (dist[i][j] && (wave_front[j] == 0) && (ind[i] != 0)){
                            wave_front[j] = iteration + 1;
                            flag = 1;
                        }
                    }
                }
            }
            iteration++;
        }
        AmountOfClusters++;
        index.push_back(iteration);
    }

    WaveClust.InitializeClusters(AmountOfClusters);

    for (i = 0; i < n; i++){
        for (j = 0; j < AmountOfClusters; j++){
            if ((wave_front[i] >= index[j]) && (wave_front[i] < index[j+1]) && (ind[i] != 0)){
                WaveClust.AddPointToCluster(j, points[i]);
            }
        }
    }
}


//Hierarchical namespace

double Hierarchical::dist(Point x1, Point x2){ return sqrt ((x1.GetX() - x2.GetX())*(x1.GetX() - x2.GetX()) + (x1.GetY() - x2.GetY())*(x1.GetY() - x2.GetY())); }
Storage Hierarchical::GetStorage(){ return storage; }

void Hierarchical::d_min(int step, int n){

    int i, j, k;
    int size = sets[step].size();
    int sizeNew;
    double d;
    int flag = 1;

    for (i = 0; i < n; i++){//пересчёт расстояний от множества до точек

        if (ind[i] == 0){//пересчитывать имеет смысл только для точек, которые ещё не учтены нигде
            d = dist (points[i], points[sets[step][0]]);
            for (j = 0; j < size; j++){
                if (d > dist(points[i], points[sets[step][j]])){ d = dist(points[i], points[sets[step][j]]); }
            }
        }

        distances[i][n + step] = d;
    }

    for (i = 0; i < step; i++){

        sizeNew = sets[i].size();
        d = dist (points[sets[i][0]], points[sets[step][0]]);
        for (j = 0; j < size; j++){
            for (k = 1; k < sizeNew; k++){
                if (d > dist(points[sets[i][k]], points[sets[step][j]])){ d = dist(points[sets[i][k]], points[sets[step][j]]); }
            }
        }//выбрано минимальное расстояние

        distances[n + i][n + step] = d;
    }

}

void Hierarchical::d_max(int step, int n){

    int i, j, k;
    int size = sets[step].size();
    int sizeNew;
    double d;
    int flag = 1;

    for (i = 0; i < n; i++){//пересчёт расстояний от множества до точек

        if (ind[i] == 0){//пересчитывать имеет смысл только для точек, которые ещё не учтены нигде
            d = dist (points[i], points[sets[step][0]]);
            for (j = 0; j < size; j++){
                if (d < dist(points[i], points[sets[step][j]])){ d = dist(points[i], points[sets[step][j]]); }
            }
        }

        distances[i][n + step] = d;
    }

    for (i = 0; i < step; i++){

        sizeNew = sets[i].size();
        d = dist (points[sets[i][0]], points[sets[step][0]]);
        for (j = 0; j < size; j++){
            for (k = 1; k < sizeNew; k++){
                if (d < dist(points[sets[i][k]], points[sets[step][j]])){ d = dist(points[sets[i][k]], points[sets[step][j]]); }
            }
        }//выбрано максимальное расстояние

        distances[n + i][n + step] = d;
    }

}

void Hierarchical::d_mid(int step, int n){}

void Hierarchical::hierarchical(Field field, int parameter){

    Point p;
    int i = 0, j = 0, k = 0;
    int IsChanged = 0;
    int n1, n2;
    double d;
    double xNew, yNew;

    points = field.GetPoints();
    int n = points.size();
    int step = n;

    ind.clear();
    ind.resize(0);
    for (i = 0; i < 2*n; i++){ ind.push_back(0); }
    distances.resize(2*n);
    sets.resize(n - 1);
    structure.clear();
    structure.resize(n - 1);
    for (i = 0; i < n; i++){ structure[i].resize(2); }

    for (i = 0; i < 2*n; i++){ distances[i].resize(2*n); }

    for (i = 0; i < n; i++){

        for (j = 0; j < n; j++){ distances[i][j] = field.distance(points[i], points[j]); }

    } //к этому моменту забита матрица расстояний

    for (k = 0; k < (n - 1); k++){//иерархическое дерево на n вершинах строится за n шагов, цикл отвечает за это

        IsChanged = 0;

        for (i = 0; i < step; i++){//тут выбираем самое маленькое расстояние из оставшихся

            for (j = (i + 1); j < step; j++){//смотрим только матрицу над диагональю -- какая-никакая, а оптимизация
                
                if ((IsChanged == 1) && (d > distances[i][j]) && (ind[i] == 0) && (ind[j] == 0)){ d = distances[i][j]; n1 = i; n2 = j; }
                if ((IsChanged == 0) && (ind[i] == 0) && (ind[j] == 0)){ d = distances[i][j]; IsChanged = 1; n1 = i; n2 = j; }

            }
        }//нашли номера "точек" с наименьшим расстоянием между ними

        ind[n1] = k + 1;// n1 < n2
        ind[n2] = k + 1;

        structure[k][0] = n1;
        structure[k][1] = n2;

        if (n1 >= n){ sets[k] = sets[n1 - n]; sets[k].insert(sets[k].end(), sets[n2 - n].begin(), sets[n2 - n].end()); }
        else if (n2 < n){sets[k].push_back(n1); sets[k].push_back(n2); }
        else { sets[k].push_back(n1); sets[k].insert(sets[k].end(), sets[n2 - n].begin(), sets[n2 - n].end()); }
//кринж

        xNew = (points[n1].GetX() + points[n2].GetX())/2.0;
        yNew = (points[n1].GetY() + points[n2].GetY())/2.0;


        p.SetX(xNew);
        p.SetY(yNew);
        points.push_back(p);//теперь надо пересчитать матрицу расстояний

        if (parameter == 1){ d_min (k, n); }//пересчитываем матрицу расстояний
        else if (parameter == 2){ d_max (k, n); }
        else if (parameter == 3){ d_mid (k, n); }


        step++;//при следующем проходе поиск осуществляем по матрице размером на 1 больше

    }

}

void Hierarchical::PrintTree(){

    ofstream file;
    ofstream sss;

    string filename;
    string pref = "name_";
    string number;
    string post = ".txt";
    stringstream s;
    string sssname = "name.txt";

    sss.open(sssname);

    int i, j;
    int n = points.size()/2;

    for (i = 1; i <= n; i++){

        s.str().resize(0);
        s.str().clear();
        s.str("");
        s << i;
        number = s.str();
        filename = pref + number + post;
        file.open(filename);

        for (j = 0; j < 2*n; j++){
            
            if (ind[j] == i) { file << points[j].GetX() << "\t" << points[j].GetY() << endl; sss << points[j].GetX() << "\t" << points[j].GetY() << endl; }

        }

        sss << endl;
        file << endl;
        file.close();
    }
}

void Hierarchical::Clusterise(int k){ //k -- количество кластеров, которые надо выделить из построенного иерархического дерева

    int n = structure.size();
    int PointsSize = points.size()/2;
    int i, j;
    int size;
    int b, number;
    int AmountOf;
    vector<int> pseudoCl; //номера множеств, являющихся кластерами
    Cluster cluster;

    PointsSize++;

    storage.ClearAllClusters();
    storage.Clear();
    storage.SetType(7);
    storage.SetSize(k);

    pseudoCl.resize(0);

    pseudoCl.push_back(structure[n - 1][0]); //
    pseudoCl.push_back(structure[n - 1][1]); //первое больше нулевого
    AmountOf = 2;

    while (AmountOf < k){//надо понять, кого дробить дальше
                         //дробить дальше надо того, кто ближе к концу (т.е. кто больше)
        b = pseudoCl[0];
        number = 0;
        for (i = 0; i < AmountOf; i++){

            if (pseudoCl[i] > b){

                b = pseudoCl[i];
                number = i;
            }
        }//нашли самого большого; теперь его надо раздробить (в смысле структуры, а не set'а)

        pseudoCl[number] = structure[b - PointsSize][0];
        pseudoCl.push_back(structure[b - PointsSize][1]);

        AmountOf++;
        
    }

    for (i = 0; i < k; i++){

        size = sets[pseudoCl[i] - PointsSize].size();

        cluster.Clear();

        for (j = 0; j < size; j++){ cluster.AddPoint(points[sets[pseudoCl[i] - PointsSize][j]]); }

        storage.AddCluster(cluster);

    }

}


//RBF namespace

vector<double> RBF::GetMS(){ return ms; }
void RBF::SetMS(vector<double> ms_){ ms = ms_; }

double RBF::Regression(vector<Point> points, vector<double> meanings, double bandwidth, Point app){

    double result = 0.0;
    double weightsum = 0.0;
    double weight;
    double exponent;
    double distance;
    double disp;
    int size;
    int i;
    int flag = 0;

    double mindist1, mindist2, mindist3;   //расстояния до трёх ближайших соседей
    int minnumber1, minnumber2, minnumber3;//номера трёх ближайших соседей

    mindist1 = d(app, points[0]);
    minnumber1 = 0;
    mindist2 = d(app, points[1]);
    minnumber2 = 1;
    mindist3 = d(app, points[2]);
    minnumber3 = 2;

    size = points.size();

    disp = bandwidth/3;//тупо от балды, правило трёх сигм и вот эта вся ерунда

    for (i = 0; i < size; i++){

        distance = d(app, points[i]);

        if (distance < bandwidth){

            exponent = (-1.0/2) * pow(distance/bandwidth, 2);//что делать с дисперсией?
            exponent /= disp;
            weight = exp(exponent);
            result += weight*meanings[i];
            weightsum += weight;
            flag = 1;

        }

        if (distance < mindist1){//на всякий случай запасаем три ближайших к точке приближения соседа

            mindist3 = mindist2;
            minnumber3 = minnumber2;
            mindist2 = mindist1;
            minnumber2 = minnumber1;
            mindist1 = distance;
            minnumber1 = i;

        }
        else if (distance < mindist2){

            mindist3 = mindist2;
            minnumber3 = minnumber2;
            mindist2 = distance;
            minnumber2 = i;

        }
        else if (distance < mindist3){

            mindist3 = distance;
            minnumber3 = i;

        }

    }

    if (!flag){ //костыль на случай отсутствия в окне точки поиска приближения соседей

        exponent = (-1.0/2) * pow(mindist1/bandwidth, 2);
        exponent /= disp;
        weight = exp(exponent);
        result += weight*meanings[minnumber1];
        weightsum += weight;

        exponent = (-1.0/2) * pow(mindist2/bandwidth, 2);
        exponent /= disp;
        weight = exp(exponent);
        result += weight*meanings[minnumber2];
        weightsum += weight;

        exponent = (-1.0/2) * pow(mindist3/bandwidth, 2);
        exponent /= disp;
        weight = exp(exponent);
        result += weight*meanings[minnumber3];
        weightsum += weight;

    }

    result /= weightsum;
//прям костыль-костыль, иначе картинка точно не получится
    if (weightsum < EPS){ result = (meanings[minnumber1] + meanings[minnumber2] + meanings[minnumber3])/3.0; }

    return result;

}

void RBF::Inter(vector<Point> points, vector<double> meanings, double bandwidth, int sampling){//sampling -- частота дискретизации

    int i, j;
    int size;
    double minX, maxX;//параметры бруса, на котором мы живём
    double minY, maxY;
    double x, y;//координаты точки, в которой будем приближать
    double value;//значение в точке, в которой приблизили
    double windowX, windowY;
    double stepX, stepY;
    Point app;
    double OneMoreValue;

    ofstream output;
    ofstream demo;
    string OutName = "interpolation_3d.txt";
    string DemoName = "smt.txt";

    output.open(OutName);
    demo.open(DemoName);

    size = points.size();
    minX = points[0].GetX(); maxX = points[0].GetX();
    minY = points[0].GetY(); maxY = points[0].GetY();

    for (i = 0; i < size; i++){

        demo << points[i].GetX() << " " << points[i].GetY() << " " << meanings[i] << endl;

        if (minX > points[i].GetX()){ minX = points[i].GetX(); }
        else if (maxX < points[i].GetX()){ maxX = points[i].GetX(); }
        if (minY > points[i].GetY()){ minY = points[i].GetY(); }
        else if (maxY < points[i].GetY()){ maxY = points[i].GetY(); }

    }

    windowX = maxX - minX;
    windowY = maxY - minY;

    stepX = windowX/double(sampling);
    stepY = windowY/double(sampling);

    x = minX; y = minY;

    for (i = 0; i < (sampling + 1); i++){

        for (j = 0; j < (sampling + 1); j++){

            app.SetX(x); app.SetY(y);
            value = Regression(points, meanings, bandwidth, app);
            app.SetX(x + stepX);
            OneMoreValue = Regression(points, meanings, bandwidth, app);
            output << x << " " << y << " " << value << endl;
            output << (x + stepX) << " " << y << " " << OneMoreValue << endl << endl << endl;

            app.SetX(x); app.SetY(y + stepY);
            OneMoreValue = Regression(points, meanings, bandwidth, app);
            output << x << " " << y << " " << value << endl;
            output << x << " " << (y + stepY) << " " << OneMoreValue << endl << endl << endl;

            x += stepX;

        }

        x = minX;
        y += stepY;

    }

    output.close();

    cout << minX << "\t" << maxX << endl;
    cout << minY << "\t" << maxY << endl;

}

double RBF::Regression1(vector<Point> points, double bandwidth, double x_){

    int i;
    int size;
    double result = 0.0;
    double weightsum = 0.0;
    double weight;
    double exponent;
    double distance;
    double disp;
    int flag = 0;

    double mindist1, mindist2; //расстояние до ближайшего соседа, расстояние до второго ближайшего соседа
    int minnumber1, minnumber2;//номер ближайшего соседа, номер второго ближайшего соседа

    mindist1 = fabs(points[0].GetX() - x_);//выбираем первичных ближайших соседей
    minnumber1 = 0;
    mindist2 = fabs(points[1].GetX() - x_);
    minnumber2 = 1;

    size = points.size();

    disp = bandwidth/3;//в качестве дисперсии берём чёт достаточно случайное, хз

    for (i = 0; i < size; i++){

        distance = fabs(points[i].GetX() - x_);

        if (distance < bandwidth){

            exponent = (-1.0/2) * pow(distance/bandwidth, 2);
            exponent /= disp;
            weight = exp(exponent);
            result += weight*points[i].GetY();
            weightsum += weight;
            flag = 1;

        }

        if (distance < mindist1){

            mindist2 = mindist1;
            minnumber2 = minnumber1;
            mindist1 = distance;
            minnumber1 = i;

        }
        else if (distance < mindist2){

            mindist2 = distance;
            minnumber2 = i;

        }

    }

    if (!flag){//костыль на случай, если в окне точки аппроксимации нет соседей (приближаем по двум ближайшим)

        exponent = (-1.0/2) * pow(mindist1/bandwidth, 2);
        exponent /= disp;
        weight = exp(exponent);
        result += weight*points[minnumber1].GetY();
        weightsum += weight;

        exponent = (-1.0/2) * pow(mindist2/bandwidth, 2);
        exponent /= disp;
        weight = exp(exponent);
        result += weight*points[minnumber2].GetY();
        weightsum += weight;

    }

    result /= weightsum;

    return result;

}

void RBF::Inter1(vector<Point> points, double bandwidth, int sampling){

    int i;
    int size;
    double step;
    double min, max;
    double window;
    double value;
    double x;

    ofstream output;
    string OutName = "interpolation.txt";

    output.open(OutName);

    size = points.size();
    min = points[0].GetX();
    max = points[0].GetX();

    for (i = 1; i < size; i++){

        if (min > points[i].GetX()){ min = points[i].GetX(); }
        if (max < points[i].GetX()){ max = points[i].GetX(); }

    }

    window = max - min;
    x = min;

    step = window/double(sampling);

    for (i = 0; i < (sampling + 1); i++){

        value = Regression1(points, bandwidth, x);
        output << x << "\t" << value << endl;
        x += step;

    }

    output.close();

}

void RBF::Generate(int amount, double center, double disp){

    int i;
    double number;
    default_random_engine gener;
    normal_distribution<double> distribution(center, disp);

    ms.clear();
    ms.resize(0);

    for (i = 0; i < amount; i++){

        number = distribution(gener);
        ms.push_back(number);

    }

}


//Exec namespace

vector<Storage> Exec::GetStorages(){ return storages; }

void Exec::AddStorage(Storage st){ storages.push_back(st); }
void Exec::SetStorages(vector<Storage> st_){ storages = st_; }

void Exec::PrintEllipses(vector<vector<vector<double>>> cov, vector<Point> centers){

    ofstream file;
    int i;
    int n = cov.size();
    double det;

    file.open("makee.txt");

    for (i = 0; i < n; i++){

        det = cov[i][0][0]*cov[i][1][1] - cov[i][1][0]*cov[i][0][1];

        file << "f" << i << "(x, y) = " << cov[i][1][1]/det << "*(x - " << centers[i].GetX() << ")**2 - ";
        file << (cov[i][1][0] + cov[i][0][1])/det << "*(x - " << centers[i].GetX() << ")*(y - ";
        file << centers[i].GetY() << ") + " << cov[i][0][0]/det << "*(y - " << centers[i].GetY();
        file << ")**2" << endl;

    }
    file.close();
}

void Exec::REplaceStorage(int i, Storage st_){ storages[i] = st_; }

void Exec::SetSaved(){

    int n = storages.size();
    int i;

    for (i = 0; i < n; i++){ storages[i].SetIsSaved(1); }

}

void Exec::Clear(){ storages.clear(); }


//DataBase_Field namespace

Field DataBase_Field::GetField(){ return field; }

void DataBase_Field::SetField(Field field_){ field = field_; }

void DataBase_Field::Save_Field(Field field){

    int number;
    int AmountOf;
    int i;

    ofstream FILEattributes;
    ofstream FILEpoints;

    string attribute_filename = "fields_attributes.txt";

	string points_filename = "\0";
    string prefix = "Field_";
    string field_number = "\0";
    string postfix = ".txt";
    stringstream n;

    vector<Point> points;

    number = field.GetID();
    AmountOf = field.GetSize();
    
    n.str().resize(0);
    n.str().clear();
    n.str("");
    n << (number);
    field_number = n.str();

    points_filename = prefix + field_number + postfix;

    FILEattributes.open(attribute_filename, ios::in | ios::app);
    FILEpoints.open(points_filename);

    points = field.GetPoints();


    FILEattributes << number << "\t" << field.GetSize() << "\t" << field.GetCenter().GetX() << "\t" << field.GetCenter().GetY() << "\t" <<
                                        field.GetEig()[0] << "\t" << field.GetEig1()[0] << "\t" << field.GetEig1()[1] << "\t" <<
                                        field.GetEig()[1] << "\t" << field.GetEig2()[0] << "\t" << field.GetEig2()[1] << "\t" <<
                                        field.GetAmountOfCl() << endl;


    for (i = 0; i < AmountOf; i++){
        FILEpoints << i << "\t" << points[i].GetX() << "\t" << points[i].GetY() << "\t" << points[i].GetColour() << endl;
    }

    FILEattributes.close();
    FILEpoints.close();

}

void DataBase_Field::UpdateAmount(){

    fstream file;
    int value;

    string filename = "fields_attributes.txt";

    file.open(filename, ios::in | ios::out);
    file >> value;
    file.seekp(0, ios::beg);
    value++;
    file << value << " -- AmountOfFieldsSaved";
    file.close();

}

int DataBase_Field::Read_Field(int number){

    Point p;
    string AttributesFileName = "fields_attributes.txt"; 
    string prefix = "Field_";
    string FieldNumber = "\0";
    string postfix = ".txt";
    string FieldFileName;
    string s;
    int AmountOf;
    int a;
    int N;
    int i, j, colour;
    int nu, FieldCard;
    double centerX, centerY;
    double x, y;
    vector <double> eigval;
    vector <double> eig1, eig2;
    stringstream n;

    ifstream AttributesFile;
    ifstream FieldFile;

    field.Clear();

    eigval.resize(2);
    eig1.resize(2);
    eig2.resize(2);

    AttributesFile.open(AttributesFileName, ios::in);
    AttributesFile >> AmountOf;

    if (number > AmountOf) return -1;

    n.str().resize(0);
    n.str().clear();
    n.str("");
    n << (number);
    FieldNumber = n.str();

    for (i = 0; i < number; i++){ getline (AttributesFile, s); }

    AttributesFile >> nu;
    AttributesFile >> N;
    AttributesFile >> centerX;
    AttributesFile >> centerY;
    AttributesFile >> eigval[0];
    AttributesFile >> eig1[0];
    AttributesFile >> eig1[1];
    AttributesFile >> eigval[1];
    AttributesFile >> eig2[0];
    AttributesFile >> eig2[1];
    AttributesFile >> a;

    p.SetX(centerX);
    p.SetY(centerY);

    field.AddCenter(p);
    field.SetEig(eigval);
    field.AddEig1(eig1);
    field.AddEig2(eig2);
    //field.SetID(AmountOf);
    field.SetID(nu);
    field.SetIsSaved(1);
    field.SetAmountOfCl(a);

    FieldFileName = prefix + FieldNumber + postfix;

    FieldFile.open(FieldFileName, ios::in);

    for (i = 0; i < N; i++){

        FieldFile >> j;
        FieldFile >> x;
        FieldFile >> y;
        FieldFile >> colour;
        getline (FieldFile, s);

        p.Set(x, y, j, colour);

        field.AddPoint(p);

    }

    FieldFile.close();
    AttributesFile.close();

    return 1;

}

int DataBase_Field::GetFieldNumber(){

    fstream file;
    int number;

    string filename = "fields_attributes.txt";

    file.open(filename, ios::in | ios::app);

    if (file.peek() == EOF){

        file.close();
        file.open(filename, ios::out);
        file << 0 << " -- AmountOfFieldsSaved" << endl;

        return 0;
    }

    file.close();
    file.open(filename, ios::in);
    file >> number;

    return number;

}


//DataBase_Cluster namespace

void DataBase_Cluster::SetCluster(Cluster cluster_){ cluster = cluster_; }
Cluster DataBase_Cluster::GetCluster(){ return cluster; }

void DataBase_Cluster::Save_Cluster(Cluster cluster_, int FNumber, int ClNumber, int EtEN, int number, int TypeOf){
                                    //кластер, который сохраняем, номер поля, номер кластеризации, номер в сквозной нумерации, номер кластера, тип кластеризации
    ofstream ClustersAttributes;//файл с атрибутами кластеров
    ofstream ClusterOut;//файл с точками конкретного кластера

    string AttName = "clusters_attributes.txt";
    string OutName;
    string dot = ".";
    string FieldNumber;
    string ClusteringNumber;
    string ClusterNumber;
    string postfix = ".txt";

    int i;
    int size = cluster_.GetPoints().size();

    stringstream f, cing, c;

    f.str().resize(0);
    cing.str().resize(0);
    c.str().resize(0);

    f.str().clear();
    cing.str().clear();
    c.str().clear();

    f.str("");
    cing.str("");
    c.str("");

    f << FNumber;
    cing << (ClNumber + 1);
    c << number;

    FieldNumber = f.str();
    ClusteringNumber = cing.str();
    ClusterNumber = c.str();

    OutName = FieldNumber + dot + ClusteringNumber + dot + ClusterNumber + postfix;

    ClustersAttributes.open(AttName, ios::in | ios::app);
    ClusterOut.open(OutName);
    ClustersAttributes << EtEN << "\t" << TypeOf << "\t" << size << "\t";

    if ((TypeOf != 2) && (TypeOf != 4)){

        ClustersAttributes << cluster_.GetCenter().GetX() << "\t" << cluster_.GetCenter().GetY() << "\t"
                           << cluster_.GetEig1()[0] << "\t" << cluster_.GetEig1()[1] << "\t"
                           << cluster_.GetEig2()[0] << "\t" << cluster_.GetEig2()[1] << "\t";

    }

    ClustersAttributes << FNumber << "\t" << ClNumber << "\t" << number << "\t" << OutName << endl;

    ClustersAttributes.close();
    
    for(i = 0; i < size; i++){ ClusterOut << cluster_.GetPoints()[i].GetNumber() << "\t" << cluster_.GetPoints()[i].GetX() << "\t" << cluster_.GetPoints()[i].GetY() << endl; }

    ClusterOut.close();

}

void DataBase_Cluster::Read_Cluster(int FNumber, int ClNumber, int number){ //номер поля, номер кластеризации (с 0), номер кластера (начиная с 1)

    ifstream ClustersAttributes;
    ifstream ClusterIn;

    string AttName = "clusters_attributes.txt";
    string InName;
    string dot = ".";
    string FieldNumber;
    string ClusteringNumber;
    string ClusterNumber;
    string postfix = ".txt";
    string s;

    Point p;

    int k;
    int i;
    int size;
    int nu;
    int TypeOf;
    int fn, cn, n;
    double centerX, centerY;
    double x_, y_;
    //double eig1, eig2;//собственные значения
    vector <double> eig1, eig2;

    stringstream f, cing, c;

    cluster.Clear();
    eig1.clear();
    eig2.clear();
    eig1.resize(0);
    eig2.resize(0);
    eig1.resize(2);
    eig2.resize(2);

    f.str().resize(0);
    cing.str().resize(0);
    c.str().resize(0);

    f.str().clear();
    cing.str().clear();
    c.str().clear();

    f.str("");
    cing.str("");
    c.str("");

    f << FNumber;
    cing << (ClNumber + 1);
    c << number;

    FieldNumber = f.str();
    ClusteringNumber = cing.str();
    ClusterNumber = c.str();

    InName = FieldNumber + dot + ClusteringNumber + dot + ClusterNumber + postfix;

    ClustersAttributes.open(AttName, ios::in | ios::app);
    ClusterIn.open(InName);
    ClustersAttributes.seekg(0, ios::beg);
    ClusterIn.seekg(0, ios::beg);

    while (!ClustersAttributes.eof()){

        ClustersAttributes >> k >> TypeOf >> size;
        if ((TypeOf != 2) && (TypeOf != 4)){ ClustersAttributes >> centerX >> centerY >> eig1[0] >> eig1[1] >> eig2[0] >> eig2[1]; }
        ClustersAttributes >> fn >> cn >> n;
        getline (ClustersAttributes, s);

        if ((fn == FNumber) && (cn == ClNumber) && (n == number)){ break; }

    }

    p.Set(centerX, centerY, n, n);
    cluster.AddCenter(p);
    cluster.AddEig1(eig1);
    cluster.AddEig2(eig2);

    for (i = 0; i < size; i++){

        ClusterIn >> nu >> x_ >> y_;
        getline (ClusterIn, s);
        p.Set(x_, y_, nu, nu);
        cluster.AddPoint(p);

    }

    ClustersAttributes.close();
    ClusterIn.close();
}

int DataBase_Cluster::UpdateAmount(){

    ifstream in;

    string InName = "clusters_attributes.txt";
    string s;

    int AmountOf = 0;

    in.open(InName, ios::in | ios::app);
    in.seekg(0, ios::beg);

    while(!in.eof()){

        getline (in, s);
        AmountOf++;

    }

    return AmountOf;
}


//DataBase_Storage namespace

DataBase_Storage::DataBase_Storage(){

    IsUpdated = 0;
    ClN = 0;
    EtEN = 0;
    AmountOf = 0;

}

void DataBase_Storage::SetStorage(Storage storage_){ storage = storage_; }
Storage DataBase_Storage::GetStorage(){ return storage; }
void DataBase_Storage::SetStorages(vector<Storage> storages_){ storages = storages_; }
vector<Storage> DataBase_Storage::GetStorages(){ return storages; }
void DataBase_Storage::SetAmountOf(int AmountOf_){ AmountOf = AmountOf_; }

void DataBase_Storage::Save_Storage(Storage storage_, int FNumber, int number){

    ofstream StorageAttributes;

    string AttName = "storage_attributes.txt";

    int type;
    int size;
    int i;

    if (storage_.GetIsSaved()){ return; }

    StorageAttributes.open(AttName, ios::in | ios::app);

    type = storage_.GetType();
    size = storage_.GetSize();

    if (!IsUpdated){ UpdateAmount(); }

    StorageAttributes << EtEN << "\t" << type << "\t";
                                              //14(?) символов
    if (type == 1){      StorageAttributes << "k-means_______ " << "\t" << size; }
    else if (type == 2){ StorageAttributes << "WaveFront_____ " << "\t" << size << "\t" << storage_.GetEps(); }
    else if (type == 3){ StorageAttributes << "WaveClustering " << "\t" << size << "\t" << storage_.GetEps(); }
    else if (type == 4){ StorageAttributes << "DB-division___ " << "\t" << size << "\t" << storage_.GetEps() << "\t" << storage_.GetF(); }
    else if (type == 5){ StorageAttributes << "DBSCAN________ " << "\t" << size << "\t" << storage_.GetEps() << "\t" << storage_.GetF(); }
    else if (type == 6){ StorageAttributes << "EM____________ " << "\t" << size; }
    else if (type == 7){ StorageAttributes << "hierarchical__ " << "\t" << size; }

    StorageAttributes << "\t" << FNumber << "\t" << (number + 1) << endl;

    StorageAttributes.close();

    for(i = 0; i < size; i++){

        dbc.Save_Cluster(storage_.GetClusters()[i], FNumber, number, ClN, (i + 1), type);
        ClN++;

    }

    EtEN++;
}

void DataBase_Storage::Read_Storage(int FNumber, int number){//номер поля, номер кластеризации (начиная с нуля)

    ifstream StorageAttributes;

    string AttName = "storage_attributes.txt";
    string s;

    int e, i;
    int type;
    int size;
    int F;
    double eps;
    int fn, cn;

    storage.Clear();

    StorageAttributes.open(AttName, ios::in);

    StorageAttributes.seekg(0, ios::beg);

    while (!StorageAttributes.eof()){

        StorageAttributes >> e >> type;
        getline(StorageAttributes, s, ' ');
        StorageAttributes >> size;
        if (type == 2){      StorageAttributes >> eps; }
        else if (type == 3){ StorageAttributes >> eps; }
        else if (type == 4){ StorageAttributes >> eps >> F; }
        else if (type == 5){ StorageAttributes >> eps >> F; }
        StorageAttributes >> fn >> cn;
        getline (StorageAttributes, s);

        if ((fn == FNumber) && (cn == (number + 1))){ break; }

    }

    storage.SetType(type);
    storage.SetSize(size);
    if ((type == 2) || (type == 3)){ storage.SetEps(eps); }
    if ((type == 3) || (type == 5)){ storage.SetEps(eps); storage.SetF(F); }

    for (i = 1; i <= size; i++){

        dbc.Read_Cluster(FNumber, number, i);//номер поля, номер кластеризации (с нуля), номер кластера (с единицы)
        storage.AddCluster(dbc.GetCluster());

    }

    StorageAttributes.close();
}

void DataBase_Storage::UpdateAmount(){

    ifstream file;
    string FileName = "storage_attributes.txt";
    string s;

    int AmountOf = 0;

    file.open(FileName, ios::in | ios::app);
    file.seekg(0, ios::beg);

    while (!file.eof()){

        getline (file, s);
        AmountOf++;

    }

    ClN = dbc.UpdateAmount();
    EtEN = AmountOf;
    IsUpdated = 1;
    
}

void DataBase_Storage::Save_Storages(vector<Storage> storages_, int FNumber){

    int i;
    int size;

    size = storages_.size();

    for (i = 0; i < size; i++){ Save_Storage(storages_[i], FNumber, i); }

}

void DataBase_Storage::Read_Storages(int FNumber){

    int i;

    storages.clear();
    storages.resize(0);

    for (i = 0; i < AmountOf; i++){

        Read_Storage(FNumber, i);
        storages.push_back(storage);
        
    }
}