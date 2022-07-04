#include<bits/stdc++.h>
#include "./christofides-algorithm/MST.h"
#include "./christofides-algorithm/Matching/Matching.h"
#include "./christofides-algorithm/Matching/Graph.h"
#include "./christofides-algorithm/Christofides.h"
#ifdef WATCH_MEM
	#include "monitor.h"
    int usedMemory;
#endif
using namespace std;
typedef pair<int, int> PI;
typedef pair<PI, double> PID;
typedef pair<double, int> P;
typedef pair<double, double> PD;
typedef pair<double, P> DDI;
const int MAX_V = 800000, MAX_M = 1500000;// V includes N and O
struct edge{
    int to, id;
    double cost;
    edge(int _to, int _id, double _cost){
        to = _to;
        id = _id;
        cost = _cost;
    }
};

//for a new or existing stop location
struct node{
    double lon, lat, dis, utility;
    double lb_dist, realcost;
    int eid, lbind;
    vector<int> routes;
    string sroutes;
    vector<P> neighbors;
    node(double _lon, double _lat, int _eid, double _dis){
        lon = _lon;
        lat = _lat;
        eid = _eid;
        dis = _dis;
        utility = 0;
        lbind = 0;
        lb_dist = 1e10;
    }
};
vector<node> NO;

double rad(double d)
{
    const double PI = 3.1415926535898;
    return d * PI / 180.0;
}

//distance between two coordinates 
double calDist(float fLati1, float fLong1, float fLati2, float fLong2)
{
    const float EARTH_RADIUS = 6378.137;
 
    double radLat1 = rad(fLati1);
    double radLat2 = rad(fLati2);
    double a = radLat1 - radLat2;
    double b = rad(fLong1) - rad(fLong2);
    double s = 2 * asin(sqrt(pow(sin(a/2),2) + cos(radLat1)*cos(radLat2)*pow(sin(b/2),2)));
    s = s * EARTH_RADIUS;
    s = (int)(s * 10000000) / 10000;
    return s/1000.0;
}

struct query{
    int num=0, nn;//nearest stop
    double nndist;//distance to the nearest stop
};

vector<edge> G1[MAX_V];//road network
vector<PD> coord;
vector<PID> edges;
int pre[MAX_V], pre_i[MAX_V], V, M, K, nN, nO, nQ;
double alpha, D, dist[MAX_V];//kilometers
query Q[MAX_V];
FILE *fp_query, *fp_transit, *fp_network;
map<string, int> r2id;
vector<string> id2r;
double service_old;
priority_queue<P> uque;
vector<int> res, busroute, busstop;
bool covered_routes[1000], T[MAX_V], validu[MAX_V], validc[MAX_V], busStop[MAX_V];
double ecost[300][300], initial_walkingDis;
int startNode;

//distance between s and t
double dijkstraCost(int s, int t){
    priority_queue<P, vector<P>, greater<P> > que;
    for (int i = 0; i < V + nN + nO; i++)
        dist[i] = DBL_MAX;
    dist[s] = 0;
    que.push(P(0, s));
    while (!que.empty()){
        P p = que.top();
        que.pop();
        int v = p.second;
        if (dist[v] < p.first)
            continue;
        if (v == t)
            return p.first;
        for (int i = 0; i < G1[v].size(); i++) {
            edge e = G1[v][i];
            if (dist[e.to] > dist[v] + e.cost){
                dist[e.to] = dist[v] + e.cost;
                que.push(P(dist[e.to], e.to));
            }
        }
    }
    return -1;
}

//build the final graph for the Christofides' algorithm
void ReadEuclideanGraph(string filename, Graph & G, vector<double> & X, vector<double> & Y, vector<double> & cost)
{
	ifstream file;
	file.open(filename.c_str());

	string s;
	getline(file, s);
	stringstream ss(s);
	int n;
	ss >> n;

    G = Graph(n);
    X.clear(); 
    Y.clear();
    cost.clear();
	for(int i = 0; i < n; i++)
	{
		getline(file, s);
		ss.str(s);
		ss.clear();
		double x, y;
		ss >> x >> y;

        X.push_back(x);
        Y.push_back(y);
	}
	for(int i = 0; i < n; i++)
	    for(int j = i+1; j < n; j++)
	        G.AddEdge(i, j);
	for(int i = 0; i < G.GetNumEdges(); i++)
	{
	    pair<int, int> p = G.GetEdge(i);
	    int u = p.first, v = p.second;
        //cost.push_back(calDist(Y[u], X[u], Y[v], X[v]));
        double c = dijkstraCost(res[u], res[v]);
        cost.push_back(c);
        ecost[u][v] = c;
    }
	file.close();
}

//nearest neighbor query called by preprocessing
void nnQ(int s){//starting from each query to N or O
    priority_queue<P, vector<P>, greater<P> > que;
    for (int i = 0; i < V + nN + nO; i++)
        dist[i] = DBL_MAX;
    dist[s] = 0;
    que.push(P(0, s));
    while (!que.empty()){
        P p = que.top();
        que.pop();
        int v = p.second;
        if (dist[v] < p.first)
            continue;
        if (v >= V + nN){
            service_old += p.first;
            Q[s].nn = v;
            Q[s].nndist = p.first;
            initial_walkingDis += p.first * Q[s].num;
            return;
        }
        if (v >= V && v < V + nN){
            NO[v - V].neighbors.push_back(P(p.first, s));
        }
        for (int i = 0; i < G1[v].size(); i++) {
            edge e = G1[v][i];
            if (dist[e.to] > dist[v] + e.cost){
                dist[e.to] = dist[v] + e.cost;
                que.push(P(dist[e.to], e.to));
            }
        }
    }
    return;
}

void add_edge(int nid, int eid, double dis0){
    PID e = edges[eid];
    int a = e.first.first, b = e.first.second;
    double dis1 = e.second - dis0;
    if(dis1<0)
        printf("Error in add_edge");
    G1[a].push_back(edge(nid, M++, dis0));
    G1[nid].push_back(edge(a, M++, dis0));
    G1[b].push_back(edge(nid, M++, dis1));
    G1[nid].push_back(edge(b, M++, dis1));
}

void takeCoord(int v, double &x, double &y){
    if (v < V){
        x = coord[v].first;
        y = coord[v].second;
    }
    else{
        x = NO[v - V].lon;
        y = NO[v - V].lat;
    }
}

void readwrite(string fs, FILE *fpw){
    FILE *fpr = fopen(fs.c_str(), "r");
    char c;
    while((c=fgetc(fpr))!=EOF)
        fputc(c, fpw);
}

//for visualization only
void plot(){
    FILE *fp_EBR = fopen("dispatch_EBR.html", "w");
    readwrite(string("plot/1.txt"), fp_EBR);
    fprintf(fp_EBR, "%f,%f", coord[0].first, coord[0].second);
    readwrite(string("plot/2.txt"), fp_EBR);
    FILE *fp_plotR = fopen("plot_route.txt", "w");
    fprintf(fp_plotR, "{'type': 'Feature','properties': {'color': '#008B8B','width': 5},'geometry': {'type': 'LineString','coordinates': [");
    fprintf(fp_EBR, "{'type': 'Feature','properties': {'color': '#008B8B','width': 5},'geometry': {'type': 'LineString','coordinates': [");
    for (int i = 0; i < busroute.size(); i++){
        int v = busroute[i];
        double x, y;
        if (v < V){
            x = coord[v].first;
            y = coord[v].second;
        }
        else{
            x = NO[v - V].lon;
            y = NO[v - V].lat;
        }
        fprintf(fp_plotR, "[%f,%f],", x, y);
        fprintf(fp_EBR, "[%f,%f],", x, y);
    }
    fprintf(fp_plotR, "]}},");
    fprintf(fp_EBR, "]}},");
    readwrite(string("plot/3.txt"), fp_EBR);
    FILE *fp_plotB = fopen("plot_bus.txt", "w");/*
    //Old bus network
    for (int i = 0; i < nO; i++){
        double x = NO[nN + i].lon, y = NO[nN + i].lat;
        fprintf(fp_plotB,"{ \"type\": \"Feature\", \"properties\": { 'description':\"%d\",'icon': 'bus'},", NO[nN+i].routes.size());
        fprintf(fp_plotB,"\"geometry\": { \"type\": \"Point\", \"coordinates\": [ %f,%f] } },\n", x, y);
    }*/
    //new bus route
    for (int i = 0; i < busstop.size(); i++){
        int v = busstop[i];
        double x, y;
        takeCoord(v, x, y);
        fprintf(fp_plotB, "{ \"type\": \"Feature\", \"properties\": { 'description':\"%d\"},", i);
        fprintf(fp_plotB,"\"geometry\": { \"type\": \"Point\", \"coordinates\": [ %f,%f] } },\n", x, y);
        fprintf(fp_EBR, "{ \"type\": \"Feature\", \"properties\": { 'description':\"%d:%d\"},", i,v);
        fprintf(fp_EBR,"\"geometry\": { \"type\": \"Point\", \"coordinates\": [ %f,%f] } },\n", x, y);
    }
    readwrite(string("plot/4.txt"), fp_EBR);

    FILE *fp_plotQ = fopen("plot_query.txt", "w");
    random_device re;
    uniform_real_distribution<double> u(-0.003, 0.003);
    for (int i = 0; i < V;i++){
        if (Q[i].num != 0){
            double cx = coord[i].first, cy = coord[i].second;
            for(int j=0;j<Q[i].num;j++){
                double x = cx+u(re);
                double y = cy+u(re);
                fprintf(fp_plotQ, "{ \"type\": \"Feature\", \"properties\": { \"id\": \"%d\", \"mag\": %f},", i, (double)Q[i].num/100); //((double)Q[i].num/100<0?0:(double)Q[i].num/10));
                fprintf(fp_plotQ,"\"geometry\": { \"type\": \"Point\", \"coordinates\": [ %f,%f, 0.0 ] } },\n", x, y);
                fprintf(fp_EBR, "{ \"type\": \"Feature\", \"properties\": { \"id\": \"%d\", \"mag\": %f},", i, 0.1); //((double)Q[i].num/100<0?0:(double)Q[i].num/10));
                fprintf(fp_EBR,"\"geometry\": { \"type\": \"Point\", \"coordinates\": [ %f,%f, 0.0 ] } },\n", x, y);
            }
        }
    }
    readwrite(string("plot/5.txt"), fp_EBR);
    /*
    FILE *fp_plotNet = fopen("plot_network.txt", "w");
    printf(" %d ", edges.size());
    for (int i = 0; i < edges.size();i++){
        int p = edges[i].first.first, q = edges[i].first.second;
        if(i%2==1)
            continue;
        fprintf(fp_plotNet, "{'type': 'Feature','properties': {'color': '#F7455D','width': 2},'geometry': {'type': 'LineString','coordinates': [");
        fprintf(fp_plotNet, "[%f,%f],[%f,%f]]}},\n", coord[p].first, coord[p].second, coord[q].first, coord[q].second);
    }*/
}

//collect the input
void input(){
    fscanf(fp_network,"%d%d", &V, &M);
    int a, b, f;
    double c, d, e;
    for (int i = 0;i < V; i++){
        fscanf(fp_network,"%lf%lf", &c, &d);
        coord.push_back(PD(c, d));
    }
    for (int i = 0; i < M; i++) {
        fscanf(fp_network,"%d%d%lf", &a, &b, &c);
        G1[a].push_back(edge(b, i, c));
        edges.push_back(PID(PI(a, b), c));
    }
    fscanf(fp_transit,"%d%d", &nN, &nO);
    for (int i = 0; i < nN; i++){
        fscanf(fp_transit, "%d%lf%lf%lf", &a, &e, &c, &d);
        NO.push_back(node(c, d, a, e));
        add_edge(V + i, a, e);
    }
    printf("%d %d\n", nN, nO);
    double maxoldU = 0;
    for (int i = 0; i < nO; i++){
        fscanf(fp_transit, "%d%lf%lf%lf", &a, &e, &c, &d);
        NO.push_back(node(c, d, a, e));
        char tmps[200];
        fscanf(fp_transit, "%s", tmps);
        string s = string(tmps)+",";
        string tmp="";
        for (int j = 0; j < s.size();j++){
            if(s[j]==','||j==s.size()-1){
                if(r2id.count(tmp)==1)
                    NO[NO.size() - 1].routes.push_back(r2id[tmp]);
                else{
                    id2r.push_back(tmp);
                    r2id[tmp] = id2r.size()-1;
                    NO[NO.size() - 1].routes.push_back(id2r.size()-1);
                }
                tmp = "";
                continue;
            }
            tmp += s[j];
        }
        add_edge(V + nN + i, a, e);
        NO[nN + i].utility = alpha * NO[nN + i].routes.size();
        if (NO[nN + i].utility > maxoldU){
            maxoldU = NO[nN + i].utility;
            startNode = V + nN + i;
        }
        NO[nN + i].sroutes = string(tmps);
        uque.push(P(NO[nN + i].utility, V + nN + i));
    }
    fscanf(fp_query,"%d", &nQ);
    for (int i = 0; i < nQ; i++){
        fscanf(fp_query,"%d", &a);
        Q[a].num++;
    }
}

//calculate the initial utility values
void preprocess(){
    for (int i = 0; i < V ; i++){
        if (Q[i].num != 0)
            nnQ(i);
    }
    for (int i = 0; i < nN; i++){
        double utmp = 0;
        for (int j = 0; j < NO[i].neighbors.size(); j++){
            P tmp = NO[i].neighbors[j];
            utmp += Q[tmp.second].num * max(0.0,(Q[tmp.second].nndist - tmp.first));//check if it >=0
            //if(Q[tmp.second].nndist<tmp.first)
            //    printf("%d %f %f Error in preprocess",tmp.second, Q[tmp.second].nndist, tmp.first);
        }
        NO[i].utility = utmp;
        uque.push(P(utmp, V + i));
    }
    //printf("%d", uque.size());
}

//real function evaluation
double realUtility(int s){//s is from N or O 
    double ret = 0;
    if (s < V + nN){
        double ub = 0;
        for (int i = 0; i < V; i++){
            ub = max(ub, Q[i].nndist);
        }
        priority_queue<P, vector<P>, greater<P> > que;
        for (int i = 0; i < V + nN + nO; i++)
            dist[i] = DBL_MAX;
        dist[s] = 0;
        que.push(P(0, s));
        while (!que.empty()){
            P p = que.top();
            que.pop();
            int v = p.second;
            if (p.first > ub)
                break;
            if (dist[v] < p.first)
                continue;
            if (Q[v].num != 0 && p.first < Q[v].nndist){
                ret += (Q[v].nndist - p.first) * Q[v].num;
            }
            for (int i = 0; i < G1[v].size(); i++) {
                edge e = G1[v][i];
                if (dist[e.to] > dist[v] + e.cost){
                    dist[e.to] = dist[v] + e.cost;
                    que.push(P(dist[e.to], e.to));
                }
            }
        }
    }
    else{
        for (int i = 0; i < NO[s - V].routes.size(); i++){
            ret += alpha * (covered_routes[NO[s - V].routes[i]] != 1);
        }
    }
    return NO[s - V].utility = ret;
}

//real function evaluation
double realCost(int s){//s is from N or O to N or O
    priority_queue<P, vector<P>, greater<P> > que;
    for (int i = 0; i < V + nN + nO; i++)
        dist[i] = DBL_MAX;
    dist[s] = 0;
    que.push(P(0, s));
    while (!que.empty()){
        P p = que.top();
        que.pop();
        int v = p.second;
        if (T[v - V]){
            return NO[s - V].realcost = max(1.0, ceil(p.first / D));
        }
        if (dist[v] < p.first)
            continue;
        for (int i = 0; i < G1[v].size(); i++) {
            edge e = G1[v][i];
            if (dist[e.to] > dist[v] + e.cost){
                dist[e.to] = dist[v] + e.cost;
                que.push(P(dist[e.to], e.to));
            }
        }
    }
    return -1;
}


//update after a stop v_i is selected
void update(int s){//s is from N or O to Q
    if (s < V + nN){
        double ub = 0;
        for (int i = 0; i < V; i++){
            ub = max(ub, Q[i].nndist);
        }
        priority_queue<P, vector<P>, greater<P> > que;
        for (int i = 0; i < V + nN + nO; i++)
            dist[i] = DBL_MAX;
        dist[s] = 0;
        que.push(P(0, s));
        while (!que.empty()){
            P p = que.top();
            que.pop();
            int v = p.second;
            if (p.first > ub)
                break;
            if (dist[v] < p.first)
                continue;
            if (Q[v].num != 0 && p.first < Q[v].nndist){
                Q[v].nn = s;
                Q[v].nndist = p.first;
            }
            for (int i = 0; i < G1[v].size(); i++) {
                edge e = G1[v][i];
                if (dist[e.to] > dist[v] + e.cost){
                    dist[e.to] = dist[v] + e.cost;
                    que.push(P(dist[e.to], e.to));
                }
            }
        }
    }
    else{
        for (int i = 0; i < NO[s - V].routes.size(); i++){
            covered_routes[NO[s - V].routes[i]] = 1;
        }
    }
}

//lazy update on U
double takeU(int s){
    if (!validu[s - V]){
        validu[s - V] = 1;
        return realUtility(s);
    }
    else
        return NO[s - V].utility;
}

//lazy update on C
double takeC(int s){
    if (!validc[s - V]){
        validc[s - V] = 1;
        return realCost(s);
    }
    else
        return NO[s - V].realcost;
}

//use lower bound for cost
double updateLb_dist(int p){
    if (validc[p - V])
        return NO[p - V].realcost;
    for (int i = NO[p - V].lbind; i < res.size();i++){
        int s = res[i];
        NO[p - V].lb_dist = min(NO[p - V].lb_dist, calDist(NO[p - V].lat, NO[p - V].lon, NO[s - V].lat, NO[s - V].lon));
    }
    NO[p - V].lbind = res.size();
    return max(1.0, NO[p - V].lb_dist / D);
}

void writeRes(){
    FILE *fp_bus = fopen("data/bus_stops.txt", "w");
    fprintf(fp_bus, "%d\n", res.size());
    for (int i = 0; i < res.size();i++){
        int s = res[i];
        fprintf(fp_bus, "%f %f\n", NO[s - V].lon, NO[s - V].lat);
    }
    fclose(fp_bus);
}

//stop selection in each iteration
void select(int start, double unused, bool vstart){ 
    // vstart: whether to choose the start node
    if (vstart){
        T[start - V] = 1;
        res.push_back(start);
        update(start);
    }
    int peeksize = 10;
    if (nN + nO < 40)
        peeksize = 1;
    while (unused > 0){
        double lb = 0, us, cs;
        vector<P> tmpNodes;
        for (int i = 0; i < peeksize; i++){
            P p = uque.top();
            uque.pop();
            int s = p.second;
            cs = takeC(s);
            if (cs > unused){
                tmpNodes.push_back(p);
                continue;
            }
            us = takeU(s);
            tmpNodes.push_back(P(us, s));
            lb = max(lb, us / cs);
            //printf("%f %f %f %d\n", us, cs, lb, s);
        }
        for (int i = 0; i < tmpNodes.size();i++)
            uque.push(tmpNodes[i]);
        priority_queue<DDI> ucque;
        while(uque.size()){
            P p = uque.top();
            if (p.first < lb)
                break;
            double tmplb = updateLb_dist(p.second);
            //printf("(%f %d)", tmplb, p.second);
            ucque.push(DDI(p.first/tmplb,p));
            uque.pop();
        }
        printf("\nucque size:%d\n", ucque.size());
        int iter = 0;
        while(ucque.size()){
            iter += 1;
            DDI ddi = ucque.top();
            ucque.pop();
            P p = ddi.second;
            int s = p.second;
            cs = takeC(s);
            //printf("%d:%f %f %d\n", iter,ucque.top().first, cs, s);
            if (!(cs < 4 && unused < 3 && iter > 50) && cs > unused){
                uque.push(p);
                continue;
            }
            us = takeU(s);
            //printf("%d:%f %f %f %f %d\n", iter,us/cs,ucque.top().first, cs, us, s);
            if (ucque.size() == 0 || us / cs >= ucque.top().first){
                //ties are broken by choosing smaller node id
                if(us/cs==ucque.top().first&& s<ucque.top().second.second){
                    ucque.push(DDI(us/cs,P(us/cs,s)));
                    break;
                }
                res.push_back(s);
                T[s - V] = 1;
                update(s);
                unused -= cs;
                printf("%f %f %f %d\n", us/cs, us, cs, s);
                break;
            }
            else
                ucque.push(DDI(us / cs, P(us, s)));
        }
        while(ucque.size()){
            uque.push(ucque.top().second);
            ucque.pop();
        }
        int resi = res[res.size() - 1];

        if (resi < V + nN){
            for (int j = 0; j < nN; j++)
                validu[j] = validc[j] = 0;
            for (int j = 0; j < nO;j++)
                validc[j + nN] = 0;
        }
        else{
            for (int j = 0; j < nN;j++)
                validc[j] = 0;
            for (int j = 0; j < nO; j++)
                validu[j + nN] = validc[j + nN] = 0;
        }
        printf("unused:%f %d\n", unused, resi);
        if(unused<0)
            break;
    }
    writeRes();
}

//find the adjacent paths called by finalpath()
void concatenatePath(int ori, int des){
    priority_queue<P, vector<P>, greater<P> > que;
    for (int i = 0; i < V + nN + nO; i++)
        dist[i] = DBL_MAX;
    double dis_cons = 0;
    fill(pre, pre + V + nN + nO, -1);
    dist[ori] = 0;
    que.push(P(0, ori));
    while (!que.empty()){
        P p = que.top();
        que.pop();
        int v = p.second;
        if (v == des){
            vector<int> path;
            for (; v != ori; v = pre[v])
                path.push_back(pre_i[v]);
            busroute.push_back(ori);
            busstop.push_back(ori);
            printf("ori%d ", ori);
            for (int i = path.size() - 1; i >= 0; i--){
                edge e = G1[v][path[i]];
                dis_cons += e.cost;
                if (dis_cons > D && e.to != des && v != ori){
                    busstop.push_back(v);
                    printf("v%d ", v);
                    dis_cons = 0;
                }
                v = e.to;
                busroute.push_back(v);
            }
            break;
        }
        if (dist[v] < p.first)
            continue;
        for (int i = 0; i < G1[v].size(); i++) {
            edge e = G1[v][i];
            if (dist[e.to] > dist[v] + e.cost){
                dist[e.to] = dist[v] + e.cost;
                pre[e.to] = v;
                pre_i[e.to] = i;
                que.push(P(dist[e.to], e.to));
            }
        }
    }
}

//find the cost between any two stops to build the final graph for the Christofides' algorithm
vector<int> res_corder;
double estimateCost(){
    double ret = 0;
    string filename = "data/bus_stops.txt";
    Graph G;
    vector<double> cost, X, Y;
    ReadEuclideanGraph(filename, G, X, Y, cost);
    pair< vector<int> , double > p = Christofides(G, cost);
    vector<int> s = p.first;
    double maxdis = 0;
    int maxind = -1;
    for (int i = 0; i < (int)s.size() - 1; i++){
        int ori = res[s[i]], des = res[s[i + 1]];
        double adjdis = dijkstraCost(ori, des); //ecost[s[i]][s[i + 1]];
        if (maxdis < adjdis){
            maxdis = adjdis;
            maxind = i + 1;
        }
    }
    res_corder.clear();
    for (int i = maxind; i < s.size(); i++)
        res_corder.push_back(res[s[i]]);
    for (int i = 1; i < maxind; i++)
        res_corder.push_back(res[s[i]]);
    for (int j = 0; j < res_corder.size() - 1; j++){
        int ori = res_corder[j], des = res_corder[j + 1];
        ret += max(1.0, ceil(dijkstraCost(ori, des) / D));
    }
    printf("Estimated cost:%f\n", 1 + ret);
    for (int j = 0; j < res_corder.size(); j++){
        int ori = res_corder[j];
    }
    return 1 + ret;
}

double minx = 190, miny = 190, maxx = -190, maxy = -190;
double expands(int v){
    double x, y;
    takeCoord(v, x, y);
    if (v >= V && v < V + nN)
        return alpha * ((x < minx) || (x > maxx) || (y < miny) || (y > maxy));
    else
        return 0;
}

//generate the final path
void finalpath(){
    for (int i = 0; i < res.size();i++){
        double x, y;
        takeCoord(res[i], x, y);
        minx = min(x, minx);
        miny = min(y, miny);
        maxx = max(x, maxx);
        maxy = max(y, maxy);
    }
    priority_queue<P> neique;
    for (int i = 0; i < res.size(); i++){
        int s = res[i];
        priority_queue<P, vector<P>, greater<P> > que;
        for (int i = 0; i < V + nN + nO; i++)
            dist[i] = DBL_MAX;
        dist[s] = 0;
        que.push(P(0, s));
        while (!que.empty()){
            P p = que.top();
            que.pop();
            int v = p.second;
            if (p.first > D)
                break;
            if (dist[v] < p.first)
                continue;
            if (v > V && T[v - V] != 1)
                neique.push(P(NO[v - V].utility + expands(v), v));
            for (int i = 0; i < G1[v].size(); i++) {
                edge e = G1[v][i];
                if (dist[e.to] > dist[v] + e.cost){
                    dist[e.to] = dist[v] + e.cost;
                    que.push(P(dist[e.to], e.to));
                }
            }
        }
    }
    double unused = K - estimateCost();
    while (unused > 0)
    {
        P p = neique.top();
        neique.pop();
        int s = p.second;
        if (T[s - V])
            continue;
        printf("%d ", s);
        double cs = takeC(s);
        res.push_back(s);
        unused -= cs;
        T[s - V] = 1;
        update(s);
    }
    writeRes();
    printf("Pruning...\n");
    while (estimateCost() > K){
        res.pop_back();
        writeRes();
    }

    string filename = "data/bus_stops.txt";
    for (int j = 0; j < res_corder.size() - 1; j++){
        int ori = res_corder[j], des = res_corder[j + 1];
        concatenatePath(ori, des);
    }
    printf("ori%d\n", res_corder[res_corder.size() - 1]);
    busstop.push_back(res_corder[res_corder.size() - 1]);

    cout << "Nodes:" << busstop.size() << endl;
    for(int i = 0; i < busstop.size(); i++){    
        //cout << busstop[i] << endl;
        busStop[busstop[i]] = 1;
    }

    FILE *fp_busroute = fopen(filename.c_str(), "w");
    for (int i = 0; i < busstop.size(); i++){
        int v = busstop[i];
        double x, y;
        if (v < V){
            x = coord[v].first;
            y = coord[v].second;
        }
        else{
            x = NO[v - V].lon;
            y = NO[v - V].lat;
        }

        fprintf(fp_busroute, "%f %f\n", x, y);
    }
    fclose(fp_busroute);
    return;
}

//compute the walking distance
double sumOfWalking(bool bus[]){
    double walkingDis = 0;
    for (int s = 0; s < V ; s++){
        if (Q[s].num != 0){
            priority_queue<P, vector<P>, greater<P> > que;
            for (int i = 0; i < V + nN + nO; i++)
                dist[i] = DBL_MAX;
            dist[s] = 0;
            que.push(P(0, s));
            while (!que.empty()){
                P p = que.top();
                que.pop();
                int v = p.second;
                if (dist[v] < p.first)
                    continue;
                if (v >= V + nN || bus[v]){
                    walkingDis += p.first * Q[s].num;
                    break;
                }
                for (int i = 0; i < G1[v].size(); i++) {
                    edge e = G1[v][i];
                    if (dist[e.to] > dist[v] + e.cost){
                        dist[e.to] = dist[v] + e.cost;
                        que.push(P(dist[e.to], e.to));
                    }
                }
            }
        }
    }
    return walkingDis;
}

void evaluate(){
    double walkingDis = sumOfWalking(busStop);
    printf("Walking Distance: %f\n",walkingDis);
    double connectivity = 0;
    for (int i = 0; i < 1000;i++)
        if(covered_routes[i]){
            connectivity += alpha;
        }
    printf("Transfer Choices: %f\n", connectivity);
    printf("Objective: %f\n", initial_walkingDis - walkingDis + connectivity);

    return;
}

void evaluate_file(string s){
    FILE *fp_alg = fopen(s.c_str(), "r");
    double buslon[300], buslat[300];
    int i = 0;
    bool alg_busStop[MAX_V], alg_routes[3000];
    memset(alg_busStop, 0, sizeof(alg_busStop));
    memset(alg_routes, 0, sizeof(alg_routes));
    while(~fscanf(fp_alg,"%lf%lf",&buslon[i],&buslat[i])){
        double mindis = 100000000;
        int nb = -1;
        for (int j = 0; j < nN + nO; j++){
            double lx = NO[j].lon, ly = NO[j].lat;
            double dis = (lx - buslon[i]) * (lx - buslon[i]) + (ly - buslat[i]) * (ly - buslat[i]);
            if (dis < mindis){
                nb = j + V;
                mindis = dis;
            }
        }
        alg_busStop[nb] = 1;
        //printf("%d %d %.9f ", busstop[i], nb, mindis);
        if (nb > V + nN){
            for (int j = 0; j < NO[nb - V].routes.size();j++){
                alg_routes[NO[nb - V].routes[j]] = 1;
                //printf("%d ",NO[nb - V].routes[j]);
            }
        }
        //cout << endl;
        i += 1;
    }
    double walkingDis = sumOfWalking(alg_busStop);
    printf("Walking Distance: %f\n", walkingDis);
    double connectivity = 0;
    for (int i = 0; i < 3000;i++)
        if(alg_routes[i]){
            connectivity += alpha;
        }
    printf("Transfer Choices: %f\n", connectivity);

    printf("Objective: %f\n", initial_walkingDis - walkingDis + connectivity);
}


int main(int argc , char * argv[]){
    string feta, ftsp;
    if (argc > 1)
        fp_network = fopen(argv[1], "r");
    else
        fp_network = fopen("data/network.txt", "r");
    if (argc > 2)
        fp_transit = fopen(argv[2], "r");
    else
        fp_transit = fopen("data/transit.txt", "r");
    if (argc > 3)
        fp_query = fopen(argv[3], "r");
    else
        fp_query = fopen("data/queries.txt", "r");
    if (argc > 4)
        alpha = stod(string(argv[4]));
    else
        alpha = 800;
    if (argc > 5)
        K = atoi(string(argv[5]).c_str());
    else
        K = 30;
    if (argc > 6)
        D = stod(string(argv[6]));
    else
        D = 2;
    if (argc > 7)
        feta = string(argv[7]);
    else
        feta = string("data/bus_stops_eta.txt");
    if (argc > 8)
        ftsp = string(argv[8]);
    else
        ftsp = string("data/bus_stops_tsp.txt");
    
    input();
    clock_t proc_s = clock();
    printf("Preprocess...\n");
    //evaluate();
    preprocess();
    printf("Select...\n");
    //startNode = uque.top().second;
    printf("startNode%d %f", startNode, NO[startNode - V].utility);
    select(startNode, ceil(2 * K / 3), 1);
    //select(uque.top().second, ceil(2 * K / 3), 0);
    printf("Final path...\n");
    finalpath();
    clock_t proc_e = clock();

    printf("\nEvaluate...\n");
    freopen("results.txt", "a", stdout);
    cout << "Nodes:" << busstop.size() << endl;
    if (argc > 7){
        printf("EBRR:\n");
        evaluate();
        printf("ETA-Pre:\n");
        evaluate_file(feta);
        printf("Vk-TSP:\n");
        evaluate_file(ftsp);
    }
    else
        evaluate();
    double proc_d = proc_e - proc_s;
    printf("%f\n", (double)proc_d / CLOCKS_PER_SEC);
    fclose(stdout);
								
    //plot();
    return 0;
}