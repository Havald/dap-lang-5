#include <ostream>
#include <iostream>
#include <fstream>
#include <ctime>
#include <algorithm>
#include <limits>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string.h>
#include <vector>
#include <assert.h>

using namespace std;

const double infinity = numeric_limits<double>::infinity();

class Node {
	friend class Graph; // Graph class maintains all nodes
public:
	size_t NumberEdgesOut() { return OutEdges.size(); };				// How many outgoing edges has this node?
	size_t NumberEdgesIn() { return InEdges.size(); };				// How many egdes point into this node?
	size_t IndexEdgeOut( size_t const WhichOne ) { return OutEdges[WhichOne]; }; 	// Tell me the index in the Graphs's memory of 
	size_t IndexEdgeIn( size_t const WhichOne ) { return InEdges[WhichOne]; };	// the edge that goes out/in from this node.
											// which is between 0 ... NumberEdgesOut/In -1
	double 	East() const { return EastKoordinate; }		// Not needed
	double 	North() const { return NorthKoordinate; }	// Not needed
private:
	// constructor (private, nodes can only be created by Graph class)
	Node(double East, double North) : EastKoordinate(East), NorthKoordinate(North) { }
	double 	EastKoordinate, NorthKoordinate; 
	vector <size_t> OutEdges;
	vector<size_t> InEdges;
};

class Edge {
	friend class Graph; // Graph class maintains edges
public:
	size_t Source() const { return MySource; } 				// returns the source node of the edge
	size_t Target() const { return MyTarget; } 				// returns the target node of the edge
	double Cost() { return DistanceIsSelected ? MyDistance : MyTravelTime; };// Depending on selected state returns time or distance
	void SetDistanceIsSelected() { DistanceIsSelected=true; };		// Select distance as cost
	void ResetDistanceIsSelected() { DistanceIsSelected=false; };		// select travel time as cost
	int 	RoadType() const { return MyRoadType; } 			// returns the road type of the edge, not needed anyhow
private:
	// private constructor, nodes can only be created by Graph class
	Edge(size_t Source, size_t Target, double TravelTime, double Distance, int RoadType ) :
		MySource(Source), MyTarget(Target), MyTravelTime(TravelTime), MyDistance(Distance), MyRoadType(RoadType) { }
	double 	MyTravelTime; 								// travel time
	double 	MyDistance;   								// spatial distance in meters
	int    	MyRoadType;       							// road type
	size_t 	MySource;
	size_t 	MyTarget;
	bool DistanceIsSelected;
};

class Graph {
	
private:
	vector<Edge> 	TheEdges;
	vector<Node> 	TheNodes;
	bool 		HasTimeInfoTooMemory;				// tells if graph has not only distance but travel time info toopublic:
	bool 		DistanceIsSelected;
public:
	size_t 	NumberOfNodes() const { return TheNodes.size(); } 	// returns the number of nodes
	size_t 	NumberOfEdges() const { return TheEdges.size(); }	// returns the number of edges
	Edge & GiveEdge(const size_t EdgeIndex) { return TheEdges[EdgeIndex]; } 
	Node & GiveNode(const size_t NodeIndex) { return TheNodes[NodeIndex]; } 	
	void 	Load(const char *filename, const bool AssumeDirected);	// loads a roadmap in TIGER format
	bool	HasTimeInfoToo() { return HasTimeInfoTooMemory; };
	void	SelectDistanceAsCost() { for( size_t i=0; i<TheEdges.size(); i++) TheEdges[i].SetDistanceIsSelected(); };
	void 	SelectTravelTimeAsCost() { 
		if (!HasTimeInfoTooMemory) throw ("No Time Info in Graph");
		for( size_t i=0; i<TheEdges.size(); i++) TheEdges[i].ResetDistanceIsSelected();
	};
	void	SelectTimeAsCost() { if (!HasTimeInfoTooMemory) throw ("No Time Info in Graph"); DistanceIsSelected=false; };
	Graph(const char *filename, bool AssumeDirected) {
		HasTimeInfoTooMemory = false;
		ifstream is(filename);
		if(is.fail()) throw( "File not found");;
		size_t NumberOfNodes, NumberOfEdges;
		is >> NumberOfNodes; // read nodes
		for(size_t i = 0; i < NumberOfNodes; ++i) {
			size_t Index;
			double East, North;
			is >> Index >> East >> North;
			if(i != Index) throw("Bad Index");
			TheNodes.push_back(Node(East,North));
		}
		is >> NumberOfEdges; // read edges
		for(size_t i = 0; i < NumberOfEdges; ++i) {
			size_t SourceId, TargetId;
			double TravelTime, Distance;
			int RoadType;
			is >> SourceId >> TargetId >> Distance >> TravelTime >> RoadType;
			if (SourceId>=TheNodes.size() || TargetId>= TheNodes.size() || SourceId<0 || TargetId<0 ) throw ("Bad Node ID");
			
						TheEdges.push_back(Edge(SourceId, TargetId, TravelTime, Distance, RoadType));
			if (!AssumeDirected) 	TheEdges.push_back(Edge(TargetId, SourceId, TravelTime, Distance, RoadType));
		}
		
		for( size_t i=0; i<TheEdges.size(); i++) {
			try {
				TheNodes[TheEdges[i].MySource].OutEdges.push_back(i);
				TheNodes[TheEdges[i].MyTarget].InEdges.push_back(i);	
			} catch( ... ) { throw "No Memory."; };
		}
			
		for( size_t i=0; i<TheEdges.size(); i++) if (TheEdges[i].MyTravelTime != TheEdges[i].MyDistance ) Graph::HasTimeInfoTooMemory=true;
	}
	
private:
	Graph(const Graph &what); 					// not implemented! forbid access.
	Graph();
	Graph &operator=(const Graph &what); 
};

class PriorityQueue {
	class HeapElement {
	friend class PriorityQueue;
	private:
		double  MyPriority; 
		size_t MyEntry;
		size_t PlaceOfIndex;
		HeapElement(double ThePriority, size_t Node ) : MyPriority(ThePriority), MyEntry(Node) { };
	};

	vector<size_t> Index;
	vector<HeapElement> TheHeapMemory;
	bool Tainted; // Once you use ExtractMinimum, heap it is tainted and you may not insert any more.
	
	void SiftUp(int i) {
		int parent = (i-1)/2;
		while(i > 0) {
			if(TheHeapMemory[Index[i]].MyPriority < TheHeapMemory[Index[parent]].MyPriority) {
				TheHeapMemory[Index[parent]].PlaceOfIndex=i;
				TheHeapMemory[Index[i]].PlaceOfIndex=parent;
				swap(Index[i],Index[parent]);
				i = parent;
				parent = (i-1)/2;
			} else return;
		}
	}
	
	void SiftDown(int i) {
		while(2*i+1 < int(Index.size())) {
			int j = 2*i+1;
			if(j+1 < int(Index.size())) if(TheHeapMemory[Index[j]].MyPriority > TheHeapMemory[Index[j+1]].MyPriority) ++j;
			if(TheHeapMemory[Index[i]].MyPriority > TheHeapMemory[Index[j]].MyPriority) {
				TheHeapMemory[Index[j]].PlaceOfIndex = i;
				TheHeapMemory[Index[i]].PlaceOfIndex = j;
				swap(Index[i],Index[j]);
				i = j;
			} else return;
		}
	}
	
	void Delete(size_t Entry) { 
		Tainted=true;
		TheHeapMemory[Index[Index.size()-1]].PlaceOfIndex = TheHeapMemory[Entry].PlaceOfIndex;
		swap(Index[Index.size()-1],Index[TheHeapMemory[Entry].PlaceOfIndex]);
		Index.pop_back();
		size_t i = TheHeapMemory[Entry].PlaceOfIndex;
		SiftUp(i); SiftDown(i);
	}
	
public:
	explicit PriorityQueue()  { Tainted=false; };	// Create Queue;

	double GetPriority(size_t Entry) { return TheHeapMemory[Entry].MyPriority; } // Gives priority of current node
	
	bool IsEmpty() const { return Index.size() == 0; } // returns true if priority queue is empty

	void Insert(double Priority) { 	// Inserts a new node. First Node gets Entry 0; secend gets 1,...
					// You are smart if you insert the node in the same sequence as in TheGraph
		assert(!Tainted);	// Once you Extracted Minimum, the Queue is tainted.
		try { Index.push_back(Index.size()); } catch (...) { throw "No Memory."; }
		try { TheHeapMemory.push_back(HeapElement(Priority,TheHeapMemory.size())); } catch (...) { throw "No Memory."; }
		TheHeapMemory[Index.size()-1].PlaceOfIndex=Index.size()-1;
		SiftUp(Index.size()-1);
	}
	
	void DecreasePriority(size_t Entry, double Priority){	// decreases the priority of Entry to new priority
		assert(Priority <= TheHeapMemory[Entry].MyPriority ); // only decreasing is implemented!
		TheHeapMemory[Entry].MyPriority=Priority;
		SiftUp(TheHeapMemory[Entry].PlaceOfIndex);
	}
	
	void ExtractMinimum(double &Priority, size_t &Entry) {	// removes the MinimumPriority from the Queue; 
									// returns its priority prio and value in value
		Priority  = TheHeapMemory[Index[0]].MyPriority;
		Entry = TheHeapMemory[Index[0]].MyEntry;
		Delete(Index[0]);					// Here the Queue is tainted.
	}
};


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Student may not edit above this line ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

double Dijkstra(size_t Start, size_t Destination, Graph &TheGraph) {
	PriorityQueue ThePath = PriorityQueue();
	double currentValue;
	size_t NodeIndex;
	for(size_t i = 0; i < TheGraph.NumberOfNodes(); i++) {
		if (i == Start) ThePath.Insert(0);
		else ThePath.Insert(infinity);
	}
	while(!ThePath.IsEmpty()) {
		ThePath.ExtractMinimum(currentValue, NodeIndex);
		Node &currentNode = TheGraph.GiveNode(NodeIndex);
		for(size_t i = 0; i < currentNode.NumberEdgesOut(); i++) {
			Edge &currentEdge = TheGraph.GiveEdge(currentNode.IndexEdgeOut(i));
			double calculateCost = currentEdge.Cost() + ThePath.GetPriority(NodeIndex);
			if(calculateCost < ThePath.GetPriority(currentEdge.Target())) {
				ThePath.DecreasePriority(currentEdge.Target(), calculateCost);
			}
		}
	}
	return ThePath.GetPriority(Destination);
}
	
void BellmanFord(size_t Start, Graph &TheGraph, vector<double> &Result, bool CheckNegativeCycles) {
	for(size_t i = 0; i < TheGraph.NumberOfNodes(); i++) {
		Result[i] = (i == Start ? 0.0 : infinity);
	}

	for(size_t i = 0; i < TheGraph.NumberOfNodes() - 1; i++) {
		for(size_t j = 0; j < TheGraph.NumberOfEdges(); j++) {
			Edge &currentEdge = TheGraph.GiveEdge(j);
			if(Result[currentEdge.Source()] != infinity 
			  && Result[currentEdge.Source()] + currentEdge.Cost() < Result[currentEdge.Target()]) {
				Result[currentEdge.Target()] = Result[currentEdge.Source()] + currentEdge.Cost();
			}
		}
	}

	if(CheckNegativeCycles) {
		for(size_t i = 0; i < TheGraph.NumberOfEdges(); i++) {
			Edge &currentEdge = TheGraph.GiveEdge(i);
			if(Result[currentEdge.Source()] != infinity 
			  && Result[currentEdge.Source()] + currentEdge.Cost() < Result[currentEdge.Target()]) {
				throw "Negative Cycle";
			}
		}
	}
}


void APSP( Graph &TheGraph, vector< vector<double> > &Result) {
	for(size_t i = 0; i < TheGraph.NumberOfNodes(); i++) {
		for(size_t j = 0; j < TheGraph.NumberOfNodes(); j++) {
			Result[i][j] = (i == j ? 0.0 : infinity);
		}
	}
	for(size_t i = 0; i < TheGraph.NumberOfEdges(); i++) {
		Edge &currentEdge = TheGraph.GiveEdge(i);
		Result[currentEdge.Source()][currentEdge.Target()] = currentEdge.Cost();
	}

	for(size_t i = 0; i < TheGraph.NumberOfNodes(); i++) {
		for(size_t j = 0; j < TheGraph.NumberOfNodes(); j++) {
			for(size_t k = 0; k < TheGraph.NumberOfNodes(); k++) {
				if(Result[j][i] + Result[i][k] < Result[j][k]) {
					Result[j][k] = Result[j][i] + Result[i][k];
				}
			}
		}
	}
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Student may not edit below this line ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


void APSPbySSSP(Graph &TheGraph, vector< vector<double> > &Result) {
	cout << "Using SSSP" << endl;
	for(size_t i=0; i<TheGraph.NumberOfNodes(); i++)
		BellmanFord(i, TheGraph, Result[i],false);
}

void PrintMatrix(vector< vector<double> > Result, size_t OutputWidth) {

	for ( size_t k=0; k<OutputWidth+1; k++ ) cout << " ";
	for ( size_t j=0; j<Result[0].size();j++) cout << setw(OutputWidth) << j << "|";
	cout << endl;
	for (size_t i=0; i<Result.size();i++) { 
		if(!(i%10)) {		

			for ( size_t k=0; k<OutputWidth+1; k++ ) cout << " ";
			for (size_t j = 0; j<Result[0].size(); j++) {
				for (size_t k = 0; k<OutputWidth; k++) cout << "-";
				cout << "+";
			}
			cout << endl;
		}
		cout  << setw(OutputWidth) << i << ":";
		for (size_t j = 0; j<Result[i].size(); j++) {
			if (Result[i][j]==infinity) cout << setw(OutputWidth) << "";
			else cout << setw(OutputWidth) << Result[i][j];
			cout << "|"; 
		}
		cout << endl;
	}
}

void PrintConnections(Graph &TheGraph ,int OutputWidth) {
	vector< vector<double> >Result; 
	try { Result.assign(TheGraph.NumberOfNodes(),vector<double>(TheGraph.NumberOfNodes(),infinity)); } 
	catch(...) { throw "No Memory."; };

	for (size_t i = 0; i < TheGraph.NumberOfNodes(); i++) Result[i][i] = 0.0;

	for(size_t e=0; e<TheGraph.NumberOfEdges(); e++ )
		Result[TheGraph.GiveEdge(e).Source()][TheGraph.GiveEdge(e).Target()]=TheGraph.GiveEdge(e).Cost();

	PrintMatrix(Result,OutputWidth);
}


int main (int argc, char* argv[]) {
	int StartIndex=0; int DestinationIndex=0;
	bool TakeTime=false; bool AssumeDirected=false; bool DirectAPSP=true; 
	bool CheckNegativeCycles=false; bool Output=false; bool DoPrintConnections=false;
	char *FileName = _strdup("");
	int OutputWidth=3;

	const char *Alert=
"FromTo\
 [ Options ] Roadmap [ StartNode [ EndNode ] ] \n\n\
You can call me with:\n\n\
  Roadmap StartNode EndNode : and I compute distance between both nodes with Dijkstra,\n\
  Roadmap StartNode         : and I compute distance from StartNode to all others (SSSP) with Bellman-Ford,\n\
  Roadmap                   : and I compute distance between all nodes (APSP) with Floyd-Warshall or\n\
                              multiple calls of SSSP, \n\
whereas\n\n\
  Roadmap   : is the roadmap in TIGER format,\n\
  StartNode : is the integer index of the start node of the path,\n\
  EndNode   : is the integer index of the end node of path.\n\n\n\
These options may occur in any order. If conlict: last option counts:\n\n\
-n        : use multiple calls of SSSP to calculate APSP (only active with APSP),\n\
-a        : calculate APSP with Bellman-Ford (default, only active with APSP),\n\
-c        : check if graph contains a negative-weight cycle (default=false, only active with SSSP),\n\
-o        : print result (default=false, active only with SSSP and APSP),\n\
-p        : print initial path costs of edges in graph (default=false),\n\
-w width  : with of output fields in distance matrix (default:3, active only width APSP or -p)\n\
-t        : measure and print time of calculation,\n\
-d        : make a directed graph from Roadmap (default=false),\n\
-u        : make a bidirectional graph from Roadmap (default, \n\
            hint: every edge in the graph is duplicated\n\
            and inserted in either direction during loading of graph).\n";

	try {
		if (argc<2) throw  Alert;
				
		bool FileNameOK=false; bool StartIndexOK=false; bool DestinationIndexOK=false;

		for(int i=1;i<argc;i++)  { // parse options
			if (!string(argv[i]).compare("-n"))  { DirectAPSP=false; continue; }
			if (!string(argv[i]).compare("-a"))  { DirectAPSP=true; continue; }
			if (!string(argv[i]).compare("-c"))  { CheckNegativeCycles=true; continue; }
			if (!string(argv[i]).compare("-o"))  { Output=true; continue; }
			if (!string(argv[i]).compare("-p"))  { DoPrintConnections=true; continue; }
			if (!string(argv[i]).compare("-t"))  { TakeTime=true; continue; }
			if (!string(argv[i]).compare("-d"))  { AssumeDirected=true; continue; }
			if (!string(argv[i]).compare("-u"))  { AssumeDirected=false;; continue; }
			if (!string(argv[i]).compare("-w"))  {
				if (++i == argc ) throw  "-w must be followed by width" ;
				try { std::istringstream(argv[i]) >> OutputWidth; }  catch(...) { throw "Malformed number for width" ; }  
				if (OutputWidth<1) throw "Width must be >0";
				continue; 
			}
			if(!FileNameOK) { FileName=argv[i]; FileNameOK=true; continue; };
			if(!StartIndexOK) { 
				try { std::istringstream(argv[i]) >> StartIndex; }  catch(...) { throw "Malformed number for Start Node" ; } 
				StartIndexOK=true; continue; 
			};
			if(!DestinationIndexOK) { 
				try { std::istringstream(argv[i]) >> DestinationIndex; }  catch(...) { throw "Malformed number for Destination Node" ; }
				DestinationIndexOK=true; continue; 
			};
		}

		if ( string(FileName).length()==0) throw "You must provide a roadmap.";
		const char *filename = FileName;
		ifstream is(filename);
		if(is.fail()) throw "File not found.";
		Graph TheGraph(filename,AssumeDirected);
		cout << "Graph has " << TheGraph.NumberOfNodes() << " nodes and " << TheGraph.NumberOfEdges() << " edges";
		if(TheGraph.HasTimeInfoToo()) cout << " with travel time info";
		cout << "." << endl;
		for ( int run=0; run<(TheGraph.HasTimeInfoToo()?2:1) ; run++) {
			double Tick;
			if (run) TheGraph.SelectTravelTimeAsCost();
			string Measurement = run?"travel time":"distance";
			if (DoPrintConnections) {
				cout << "Initial connections (" << Measurement << "), vertical <from> Horizontal <to>:" << endl;
				PrintConnections(TheGraph,OutputWidth);
			}
			
			if (StartIndexOK) {
				if ( StartIndex<0 || StartIndex >= int(TheGraph.NumberOfNodes()) )  throw "Start node is not in graph!";
	
				if (DestinationIndexOK) { // Must be Dijkstra

					if ( DestinationIndex<0 || DestinationIndex >= int(TheGraph.NumberOfNodes()) )  throw "Destination node is not in graph!";
					if(!run) cout << "Compute distance from " << StartIndex << " to node " << DestinationIndex << "." << endl;
					Tick = double(clock());
					double Result = Dijkstra(StartIndex, DestinationIndex, TheGraph);
					Tick=(clock()-Tick)/CLOCKS_PER_SEC;
					cout 	<< "Shortest " << Measurement << " is: ";
					if (Result==infinity) cout << "no path found."; else cout << Result;   
					cout << endl;

				} else { // Do Source-Shortest-Path-Problem (SSSP) with Bellman-Ford.

					cout 	<< "Compute shortest " << Measurement 
						<< " from node " << StartIndex << " to all other nodes." << endl;
					Tick = double(clock());
					vector<double> Result; 
					try { Result.assign(TheGraph.NumberOfNodes(),false); } catch(...) { throw "No Memory."; };
					BellmanFord(StartIndex, TheGraph, Result,CheckNegativeCycles);
					Tick=(clock()-Tick)/CLOCKS_PER_SEC;
					if(Output) {
						for (size_t i = 0; i<TheGraph.NumberOfNodes(); i++)  {
							cout << setw(OutputWidth) << i << ": " << setw(OutputWidth); 
							if (Result[i]==infinity) cout << "No path found."; else cout << Result[i]; 
							cout << endl;
						}
					}
				}
			} else {
				cout << "Compute shortest " << Measurement << " between all nodes." << endl;
				Tick = double(clock());
				vector< vector<double> >Result; 
				try { Result.assign(TheGraph.NumberOfNodes(),vector<double>(TheGraph.NumberOfNodes())); } 
				catch(...) { throw "No Memory."; };

				if(DirectAPSP) APSP(TheGraph, Result); else APSPbySSSP(TheGraph, Result);

				Tick=(clock()-Tick)/CLOCKS_PER_SEC;
				if (Output) {
					cout << "vertical <from> Horizontal <to>:" << endl;
					PrintMatrix(Result,OutputWidth);
				}
			}
			if (TakeTime) cout << "It took " << Tick << " Seconds to compute." << endl;
		}
		return 0;

	} catch (const char *Reason) { cerr << Reason << endl; }
	return 1;
}

