#include "headers.hpp"


int main (void){

    Interface interface;
    Cloud cloud_;
    Storage st_;
    Cluster cl_;
    Point c;
    string s, hstring;
    ifstream Cmd;
    ifstream Help;
    ofstream file;
    ofstream log_interface, log_controller;
    ofstream info;
    vector<Point> points;
    int command_code;
	double centerX, centerY, dispX, dispY;
	double x, y, angle;
	int n, m, ii, colour, i, k, number;
    int b;
	int mode;
    int p = 0;
    int parameter;
    double threshold;
    double result;

    int AmountOf;
    
    cout << "This program can clasterize a field of points." << endl;
	cout << "Enter 1 if you want to read commands from the file Do_magic." << endl;
    cout << "Enter 2 if you want to enter commands by yourself." << endl;

    cin >> mode;
 
    log_interface.open("log_interface.txt");
    log_controller.open("log_controller.txt"); 
    if (mode){ Cmd.open("Do_magic.txt"); }

    while(1){
    
        if (mode == 2){ cin >> s; }
        
        if (mode == 1){ getline (Cmd, s);  }
        
        command_code = interface.Start(s);

        //вывести help.txt
        if (command_code == 0){

            Help.open("help.txt");
            while (getline (Help, hstring)){ cout << hstring << endl; }
            cout << endl;
            Help.close();
        
            log_controller << "The file help.txt was printed." << endl;
        }

        //создать облако
        else if (command_code == 1){

            log_interface << "The command Create_a_cloud started." << endl;

            if (mode == 2){

                cin >> n;
                cin >> centerX;
                cin >> centerY;
                cin >> dispX;
                cin >> dispY;
                cin >> colour;
            
            }

            if (mode == 1){

                Cmd >> n;
                Cmd >> centerX;
                Cmd >> centerY;
                Cmd >> dispX;
                Cmd >> dispY;
                Cmd >> colour;
                getline (Cmd, s);

            }

            points = interface.controller.Cloud_Gener(n,centerX, centerY, dispX, dispY, colour);
            cloud_.AddPoints(points);
            interface.controller.buffer.FindFactors(points);
            c.SetX(interface.controller.buffer.GetCenter()[0]);
            c.SetY(interface.controller.buffer.GetCenter()[1]);
            cloud_.AddCenter(c);

            if (!interface.controller.field.GetCondition()){
                
                interface.controller.field.AddCloud(cloud_);

                if (interface.controller.field.GetIsSaved()){

                    interface.controller.field.SetIsSaved(0);
                    interface.controller.field.SetID(interface.controller.field.GetID() + 1);

                }
            
                p++;

                log_controller << "A cloud number " << p << " was generated." << endl;
                log_controller << "The cloud was added to the field." << endl;

            }

            if (interface.controller.field.GetCondition()){ log_controller << "error" << endl; } //написать тут про то, что кластеризация уже начата и нечего добавлять чёт ещё
        }

        //раньше печатало облака, теперь ничего не делает (да оно и раньше не особо что-то делало, если честно...)
        else if (command_code == 2){}
        //    cin >> colour;
        //    i = interface.controller.Cloud_Find(colour);
        
        //печатает поле
        else if (command_code == 3){ interface.controller.Field_Print(); }

        //k-means
        else if (command_code == 4){

            log_interface << "K_means algorithm started." << endl;

            if (mode == 2){
                cout << "???" << endl;
                cin >> k;
            }
            if (mode == 1){
                Cmd >> k;
                getline (Cmd, s);
            }

            log_controller << "K-means clasterization done." << endl;

            //interface.controller.Field_Print();

            interface.controller.exec.k_means.K_means_(k, interface.controller.field);
            interface.controller.exec.AddStorage(interface.controller.exec.k_means.GetStorage());
            interface.controller.field.SetCondition(1);

            b = interface.controller.exec.GetStorages().size();
            b--;
            m = interface.controller.exec.GetStorages()[b].GetSize();
            st_ = interface.controller.exec.GetStorages()[b];
            st_.ClearAllClusters();
            st_.Clear();

            for (i = 0; i < m; i++){

                cl_ = interface.controller.exec.GetStorages()[b].GetClusters()[i];
                interface.controller.buffer.FindFactors(cl_.GetPoints());
                c.SetX(interface.controller.buffer.GetCenter()[0]);
                c.SetY(interface.controller.buffer.GetCenter()[1]);
                cl_.AddCenter(c);
                cl_.AddEig1(interface.controller.buffer.GetEig1());
                cl_.AddEig2(interface.controller.buffer.GetEig2());
                st_.AddCluster(cl_);

            }

            interface.controller.exec.REplaceStorage(b, st_);
        }

        //распечатывает кластеры
        else if (command_code == 5){

            log_interface << "Printing clusters started." << endl;

            if (mode){ Cmd >> k; getline (Cmd, s); }
            interface.controller.Print_to_File(k - 1);

            log_controller << "Clusters were printed into files." << endl;
            
        }

        //волновой алгоритм
        else if (command_code == 6){

            log_interface << "Wave algorithm started." << endl;

            if (mode == 0){ 
                cout << "???" << endl;
                cin >> threshold;
            }
            if (mode == 1){
                Cmd >> threshold;
                getline (Cmd, s);
            }

            interface.controller.exec.wave.Wave_(threshold, interface.controller.field);
            interface.controller.exec.AddStorage(interface.controller.exec.wave.GetWaveFront()); //вывод фронта волны
            //interface.controller.exec.AddStorage(interface.controller.exec.wave.GetStorage()); //вывод именно кластеризации
            interface.controller.field.SetCondition(1);
            
            log_controller << "Wave clasterization done." << endl;

            b = interface.controller.exec.GetStorages().size();
            b--;
            m = interface.controller.exec.GetStorages()[b].GetSize();
            st_ = interface.controller.exec.GetStorages()[b];
            st_.ClearAllClusters();
            st_.Clear();

            for (i = 0; i < m; i++){

                cl_ = interface.controller.exec.GetStorages()[b].GetClusters()[i];
                interface.controller.buffer.FindFactors(cl_.GetPoints());
                c.SetX(interface.controller.buffer.GetCenter()[0]);
                c.SetY(interface.controller.buffer.GetCenter()[1]);
                cl_.AddCenter(c);
                cl_.AddEig1(interface.controller.buffer.GetEig1());
                cl_.AddEig2(interface.controller.buffer.GetEig2());
                st_.AddCluster(cl_);

            }

            interface.controller.exec.REplaceStorage(b, st_);
        }

        //DBSCAN
        else if (command_code == 7){

            log_interface << "DBScan algorithm started." << endl;

            if (mode == 1){

                Cmd >> k;
                getline (Cmd, s);
                Cmd >> threshold;
                getline (Cmd, s);

            }

            interface.controller.exec.db.DBSCAN_(k, threshold, interface.controller.field);
            interface.controller.exec.AddStorage(interface.controller.exec.db.GetStorage());
            //interface.controller.exec.AddStorage(interface.controller.exec.db.GetWaveClust());
            interface.controller.field.SetCondition(1);

            log_controller << "DBScan clasterization done." << endl;

            b = interface.controller.exec.GetStorages().size();
            b--;
            m = interface.controller.exec.GetStorages()[b].GetSize();
            st_ = interface.controller.exec.GetStorages()[b];
            st_.ClearAllClusters();
            st_.Clear();

            for (i = 0; i < m; i++){

                cl_ = interface.controller.exec.GetStorages()[b].GetClusters()[i];
                interface.controller.buffer.FindFactors(cl_.GetPoints());
                c.SetX(interface.controller.buffer.GetCenter()[0]);
                c.SetY(interface.controller.buffer.GetCenter()[1]);
                cl_.AddCenter(c);
                cl_.AddEig1(interface.controller.buffer.GetEig1());
                cl_.AddEig2(interface.controller.buffer.GetEig2());
                st_.AddCluster(cl_);

            }

            interface.controller.exec.REplaceStorage(b, st_);
        }

        //аффинные преобразования
        else if (command_code == 8){

            if (mode == 1){

                Cmd >> number;
                getline (Cmd, s);
                Cmd >> angle;
                getline (Cmd, s);

            }

            log_interface << "The cloud number " << number << " was chosen to rotate." << endl;

            cloud_ = interface.controller.field.GetCloudByNumber(number);
            interface.controller.buffer.AddCloud(cloud_);
            interface.controller.buffer.Rotate(angle);
            cloud_ = interface.controller.buffer.GetCloud();
            interface.controller.field.ReplaceCloudByNumber(number, cloud_);
            interface.controller.field.UpdatePoints();

            if (interface.controller.field.GetIsSaved()){

                interface.controller.field.SetIsSaved(0);
                interface.controller.field.SetID(interface.controller.field.GetID() + 1);
            
            }

            log_controller << "The cloud number " << number << " was rotated." << endl;
        }

        //завершение программы
        else if (command_code == 9){

            log_interface << "Stop program." << endl;

            if (mode == 1){ Cmd.close(); }
            log_controller << "Program completed." << endl;
            log_interface.close();
            log_controller.close(); 
            return 0;

        }

        //минимальное покрывающее дерево
        else if (command_code == 10){

            log_interface << "Minimal spanning tree algorithm started." << endl;
            interface.controller.exec.sp_tree.SpanningTree_(interface.controller.field);
            interface.controller.field.SetTree(interface.controller.exec.sp_tree.GetIncidence());
            log_controller << "Minimal spanning tree algorithm done." << endl;

        }

        //печатает минимальное покрывающее дерево
        else if (command_code == 11){

            log_controller << "Minimal spanning tree algorithm printed." << endl;
            interface.controller.Print_Tree();

        }

        //em алгоритм
        else if (command_code == 12){

            if (mode == 0){

                cout << "???" << endl;
                cin >> k;

            }

            if (mode == 1){

                Cmd >> k;
                getline (Cmd, s);

            }

            log_interface << "EM algorithm started." << endl;
            interface.controller.exec.em.EM_(k, interface.controller.field);
            interface.controller.exec.AddStorage(interface.controller.exec.em.GetStorage());
            interface.controller.field.SetCondition(1);

            log_controller << "EM clasterization done." << endl;

            b = interface.controller.exec.GetStorages().size();
            b--;
            m = interface.controller.exec.GetStorages()[b].GetSize();
            st_ = interface.controller.exec.GetStorages()[b];
            st_.ClearAllClusters();
            st_.Clear();

            for (i = 0; i < m; i++){

                cl_ = interface.controller.exec.GetStorages()[b].GetClusters()[i];
                interface.controller.buffer.FindFactors(cl_.GetPoints());
                c.SetX(interface.controller.buffer.GetCenter()[0]);
                c.SetY(interface.controller.buffer.GetCenter()[1]);
                cl_.AddCenter(c);
                cl_.AddEig1(interface.controller.buffer.GetEig1());
                cl_.AddEig2(interface.controller.buffer.GetEig2());
                st_.AddCluster(cl_);

            }

            interface.controller.exec.REplaceStorage(b, st_);

        }

        //печатает области целостности (?)
        else if (command_code == 13){

            log_controller << "Spheres in EM algorithm printed." << endl;
            interface.controller.exec.PrintEllipses(interface.controller.exec.em.GetCov(), interface.controller.exec.em.GetCentroids());

        }

        //печатает факторы
        else if (command_code == 14){

            file.open("ff.txt");
            interface.controller.buffer.FindFactors(interface.controller.field.GetPoints());
            file << "set arrow from " << interface.controller.buffer.GetCenter()[0] << "," << interface.controller.buffer.GetCenter()[1] << " to ";
            file << interface.controller.buffer.GetEig1()[0] + interface.controller.buffer.GetCenter()[0] << "," << interface.controller.buffer.GetEig1()[1] + interface.controller.buffer.GetCenter()[1] << " lw 2" << endl;
            file << "plot (0,0) lc 0" << endl;
            file << "set arrow from " << interface.controller.buffer.GetCenter()[0] << "," << interface.controller.buffer.GetCenter()[1] << " to ";
            file << interface.controller.buffer.GetEig2()[0] + interface.controller.buffer.GetCenter()[0] << "," << interface.controller.buffer.GetEig2()[1] + interface.controller.buffer.GetCenter()[1] << " lw 2" << endl;
            file << "plot (0,0) lc 0" << endl;
            file << endl;
            log_controller << "Factors printed.";

        }

        //info
        else if (command_code == 15){

            info.open("Info.txt");
            info << "Clusterizations and their parameters." << endl << endl;

            m = interface.controller.exec.GetStorages().size();

            for (i = 0; i < m; i++){

                ii = interface.controller.exec.GetStorages()[i].GetType();

                if (ii == 1){

                    info << "K_means algorithm." << endl;
                    info << "Amount of clusters: " << interface.controller.exec.GetStorages()[i].GetSize() << endl << endl;

                }

                else if (ii == 2){

                    info << "Wave_front" << endl;
                    info << "Threshold: " << interface.controller.exec.GetStorages()[i].GetEps() << endl;
                    info << "Iterations: " << interface.controller.exec.GetStorages()[i].GetSize() << endl << endl;

                }

                else if (ii == 3){

                    info << "Wave algorithm." << endl;
                    info << "Threshold: " << interface.controller.exec.GetStorages()[i].GetEps() << endl;
                    info << "Amount of clusters: " << interface.controller.exec.GetStorages()[i].GetSize() << endl << endl;
                    
                }

                else if (ii == 4){
//вот следующие две штуки стоит както разделить с точки зрения текстовых именований
                    info << "DBSCAN algorithm." << endl;//это разбиение на core-периферию-шум
                    info << "Threshold: " << interface.controller.exec.GetStorages()[i].GetEps() << endl;
                    info << "Amount of clusters: " << interface.controller.exec.GetStorages()[i].GetSize() << endl;
                    info << "Neighbours: " << interface.controller.exec.GetStorages()[i].GetF() << endl << endl;
                    
                }

                else if (ii == 5){

                    info << "DBSCAN algorithm." << endl;//это кластеризация волной по всему кроме шума
                    info << "Threshold: " << interface.controller.exec.GetStorages()[i].GetEps() << endl;
                    info << "Amount of clusters: " << interface.controller.exec.GetStorages()[i].GetSize() << endl;
                    info << "Neighbours: " << interface.controller.exec.GetStorages()[i].GetF() << endl << endl;
                    
                }

                else if (ii == 6){

                    info << "EM algorithm." << endl;
                    info << "Amount of clusters: " << interface.controller.exec.GetStorages()[i].GetSize() << endl << endl;
                    
                }
            }

            info.close();
            cout << "Information about clasterizations is in file info.txt" << endl;
            log_interface << "Information about clasterizations was printed into info.txt" << endl;
        }

        //занести поле в картотеку
        else if (command_code == 16){

            if (interface.controller.field.GetID() < 0){ interface.controller.field.SetID(interface.controller.dbf.GetFieldNumber() + 1); }
            if (interface.controller.field.GetIsSaved()){ interface.controller.field.SetID(interface.controller.dbf.GetFieldNumber()); }// + 1);                   
            interface.controller.buffer.FindFactors(interface.controller.field.GetPoints());
            c.SetX(interface.controller.buffer.GetCenter()[0]);
            c.SetY(interface.controller.buffer.GetCenter()[1]);
            interface.controller.field.AddCenter(c);
            interface.controller.field.SetEig(interface.controller.buffer.GetEig());
            interface.controller.field.AddEig1(interface.controller.buffer.GetEig1());
            interface.controller.field.AddEig2(interface.controller.buffer.GetEig2());
            interface.controller.field.SetAmountOfCl(interface.controller.exec.GetStorages().size()); //с этим параметром ещё стоит пошаманить чуть-чуть, чтоб можно было сквозную нумерацию делать без предварительной загрузки уже сохранённого
            if (interface.controller.field.GetIsSaved()){ interface.controller.field.SetID(interface.controller.field.GetID() + 1); }
            interface.controller.dbf.Save_Field(interface.controller.field);
            interface.controller.dbf.UpdateAmount();
            interface.controller.field.SetIsSaved(1);
            log_controller << "Field saved." << endl;

            //int correct = 1;
            //cout << correct++ << "\t";

        }

        //очистить поле
        else if (command_code == 17){

            interface.controller.field.Clear();
            interface.controller.exec.Clear(); //потому что на новом поле будут уже новые кластеризации

            if (interface.controller.field.GetIsSaved()){//чтобы при очистке поля номер увеличивался на единицу 
                                                         //только если очищенное было сохранено
                interface.controller.field.SetID(interface.controller.field.GetID() + 1);
                interface.controller.field.SetIsSaved(0);
                interface.controller.field.SetAmountOfCl(0);
        
            }

            log_controller << "Field cleared." << endl;

        }

        //считать поле из картотеки
        else if (command_code == 18){

            if (mode){

                Cmd >> number;
                log_interface << "Field number " << number << " have to be read." << endl;
                getline (Cmd, s);

            }

            if (interface.controller.dbf.Read_Field(number) < 0) {

                log_controller << "There is no field under the number " << number << endl;
                break;

            }    

            else{
                
                interface.controller.SetField(interface.controller.dbf.GetField()); 
                //interface.controller.field.SetID(interface.controller.field.GetID() + 1); 
                log_controller << "Field number " << number  << " was read." << endl;

            }          

        }

        //сохранение "FindCluster'a"
        else if (command_code == 19){

            if (interface.controller.field.GetIsSaved()){
                
                interface.controller.dbs.Save_Storages(interface.controller.exec.GetStorages(), interface.controller.field.GetID());
                interface.controller.exec.SetSaved();
                log_controller << "storage_saved" << endl;
            
            }

            else { log_controller << "а как сохранять, если после не" << endl; }

        }

        //очистка "FindCluster'a"
        else if (command_code == 20){

            interface.controller.exec.Clear();
            log_controller << "Exec почистили" << endl;

        }

        //чтение "FindCluster'a"
        else if (command_code == 21){

            //Cmd >> number;
            //getline (Cmd, s);

            interface.controller.dbs.SetAmountOf(interface.controller.field.GetAmountOfCl());
            interface.controller.dbs.Read_Storages(interface.controller.field.GetID());
            interface.controller.exec.SetStorages(interface.controller.dbs.GetStorages());

        }

        //Иерархический алгоритм
        else if (command_code == 22){

            Cmd >> parameter;
            getline(Cmd, s);
            Cmd >> k;
            getline(Cmd, s);

            //cout << "a" << endl;
            interface.controller.exec.hc.hierarchical(interface.controller.field, parameter);
            //cout << "оно по крайней мере заканчиваетсz" << endl;
            interface.controller.exec.hc.PrintTree();
            //cout << "printed" << endl;
            interface.controller.exec.hc.Clusterise(k);
            //cout << "clusterised" << endl;
            interface.controller.exec.AddStorage(interface.controller.exec.hc.GetStorage());
            //cout << "Ehf" << endl;

            log_controller << "hier done" << endl;

        }

        //считать конкретное хранилище
        else if (command_code == 23){

            Cmd >> number;
            getline(Cmd, s);

            interface.controller.dbs.Read_Storage(interface.controller.field.GetID(), number - 1);
            interface.controller.exec.AddStorage(interface.controller.dbs.GetStorage());

        }

        //непараметрическая регрессия
        else if (command_code == 24){

            Cmd >> threshold;
            getline(Cmd, s);
            Cmd >> x >> y;
            getline(Cmd, s);
            c.SetX(x); c.SetY(y);

            points = interface.controller.field.GetPoints();
            interface.controller.exec.rbf.Generate(points.size(), 10, 2);
            result = interface.controller.exec.rbf.Regression(points, interface.controller.exec.rbf.GetMS(), threshold, c);

            cout << result << endl;

        }

        //регрессия: нарисовать картинку (хз насколько адекватно оно в 3d смотреться будет)
        else if(command_code == 25){

            Cmd >> threshold;
            getline(Cmd, s);
            Cmd >> parameter;
            getline(Cmd, s);

            points = interface.controller.field.GetPoints();
            interface.controller.exec.rbf.Generate(points.size(), 10, 2);
            interface.controller.exec.rbf.Inter(points, interface.controller.exec.rbf.GetMS(), threshold, parameter);

        }

        //одномерная непараметрическая регрессия
        else if(command_code == 26){

            Cmd >> threshold;
            getline(Cmd, s);
            Cmd >> x;
            getline(Cmd, s);

            points = interface.controller.field.GetPoints();
            result = interface.controller.exec.rbf.Regression1(points, threshold, x);

            cout << result << endl;
            
        }

        //нарисовать картинку под одномерную регрессию
        else if(command_code == 27){

            Cmd >> threshold;
            getline(Cmd, s);
            Cmd >> parameter; //частота дискретизации (сколько раз приближать на имеющемся отрезке)
            getline(Cmd, s);

            points = interface.controller.field.GetPoints();
            interface.controller.exec.rbf.Inter1(points, threshold, parameter);
            //cout << "11111" << endl;

        }

        else {

            if (mode == 1){ Cmd.close(); }
            cout << "Unknown command. Check log_interface.txt to see more." << endl;
            log_interface << "Unknown command. Check help.txt to see the list of available commands and rewrite Do_magic.txt." << endl;
            log_interface.close();
            log_controller.close(); 
            return 0;

        }
    }

    log_interface.close();
    log_controller.close(); 
    if (mode == 1){ Cmd.close(); }
}