#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <utility>
#include <queue>
#include <functional>
#include "slist.h"

using namespace std;

class Airport;

double distanceEarth(Airport* air1, Airport* air2);
double distanceEarth(double lat1d, double lon1d, double lat2d, double lon2d);
Slist* distanceList(Slist* list, Airport* ref, double distance);
Airport* farthest(Slist* list, Airport* ref);

void heapSort(Slist* list, int n, Airport* reference);
void heapify(Slist* list, Airport* reference, int n, int i);

void simpleSortTotal(Slist* list, Airport* reference, int low, int high);
void merge(Slist *list, Airport* reference, int low, int m, int high);

void printCode(Slist *list, Airport *reference);

class Airport: public CustomNode
{
public:
	Airport(string c, double lo, double la) : code(std::move(c)), longitude(lo), latitude(la) {}
	Airport() = default;

	string    code;
	double   longitude{};
	double   latitude{};
	string toString() override {
		return code + " long: " + to_string(longitude) + " lat: " + to_string(latitude);
	}
	string getCode() {
		return code;
	}
	double getLong() {
		return longitude;
	}
	double getLat() {
		return latitude;
	}
	void setCode(string c) {
		code = c;
	}
	void setLong(double l) {longitude = l;}
	void setLat(double l) {latitude = l;}
	void swapData(CustomNode* node2) {
		Airport* node = dynamic_cast<Airport*>(node2);
		double tempLong = longitude;
		double tempLat = latitude;
		string tempCode = code;
		longitude = node->getLong();
		latitude = node->getLat();
		setCode(node->getCode());
		node->setLat(tempLat);
		node->setCode(tempCode);
		node->setLong(tempLong);
	}
	CustomNode* copyData() override {
		auto *air = new Airport;
		air->setCode(code);
		air->setLat(latitude);
		air->setLong(longitude);
		return air;
	};
	CustomNode* copyDataAndPtr() override {
		auto *air = copyData();
		air->ptr = ptr;
		return air;
	}
	int compare(Airport* other, Airport* reference) {
		double thisD = distanceEarth(reference, this);
		double otherD = distanceEarth(reference, other);
		if(thisD == otherD) return 0;
		else if(thisD < otherD) return -1;
		else return 1;
	}
	string toString(CustomNode* ref) override {
		return toString() + " distance: " + to_string(distanceEarth(this, reinterpret_cast<Airport *>(ref)));
	}
	bool operator ==(Airport* air) {
		return air->getLong() == getLong() && air->getLat() == getLat() && getCode() == air->getCode();
	}
};

int main()
{
	ifstream infile;
	int i=0;
	char cNum[10] ;
	//Airport* airportArr[13500];			// Replace array with Linked List
	auto * airportArr = new Slist();
	Airport* aus;
	int   airportCount;
	//Airport* a[13500];
	infile.open ("C:/Users/Alex/ClionProjects/LASAlle/Airport/USAirportCodes.csv", ifstream::in);
	infile.seekg(30);
	if (infile.is_open())
    {
	    int   c=0;
	    while (infile.good())
        {
	        auto *air = new Airport();
	        char readCode[5];
	        infile.getline(readCode, 256, ',');
	        string codeStr(readCode);
	        air->code = codeStr;
	        infile.getline(cNum, 256, ',');
	        air->longitude = atof(cNum);
	        infile.getline(cNum, 256, '\n');
	        air->latitude = atof(cNum);
	        if (!(c % 1000)) {
		        cout << air->toString() <<  endl;
	        }
	        airportArr->add(air);
	        if(air->getCode() == "AUS") {
	        	aus = dynamic_cast<Airport *>(air->copyData());
	        }
	        /*
			if (!(c % 1000))
			{
				cout << airportArr[c]->code << " long: " << airportArr[c]->longitude << " lat: " << airportArr[c]->latitude <<  endl;
				cout << airportArr[c+1]->code << endl; //" long: " << airportArr[c+1]->longitude << " lat: " << airportArr[c+1]->latitude <<  endl;
			}
			*/
	        i++ ;
	        c++;
        }
	    //airportArr->toString(aus);
	    cout << "You may experience some delays while the list is sorted.\n";
	    heapSort(airportArr, airportArr->size(), aus);
	    Slist* hundred = distanceList(airportArr, aus, 100);
	    hundred->remove(0); //0 Will be AUS
	    Airport* far = dynamic_cast<Airport *>(airportArr->get(airportArr->size()-2));
	    //simpleSortTotal(airportArr, aus, 0, airportArr->size() - 1);
	    cout << "Sorting done" << endl;
	    //cout << airportArr->toString(aus) << endl;
	    airportCount = c-1;
	    infile.close();
	    for (int k=0; k < airportCount; k++) {
		    if (!(k % 1000)) {
			    cout << airportArr->get(k)->toString() << endl;
			    cout << airportArr->get(k + 1)->toString() << endl;
			    cout << "Distance between " << dynamic_cast<Airport*>(airportArr->get(k))->getCode() << " and "
			         << dynamic_cast<Airport*>(airportArr->get(k + 1))->getCode() << " is "
			         << distanceEarth(reinterpret_cast<Airport *>(airportArr->get(k)),
			                          reinterpret_cast<Airport *>(airportArr->get(k + 1))) << endl;
		    }
	    }
	    cout << "\nFarthest airport from AUS:\n";
	    cout << far->toString(aus) << endl;
	    cout << "\nAirports within 100 miles from AUS:\n";
	    hundred->toString(aus);
	    delete hundred;
	    delete aus;
	    delete far;
	    delete airportArr;
    } else {
	    cout << "Error opening file";
    }
}

#define pi 3.14159265358979323846

#define earthRadiusKm 6371.0

// This function converts decimal degrees to radians
double deg2rad(double deg) {
	return (deg * pi / 180);
}

//  This function converts radians to decimal degrees
double rad2deg(double rad) {
	return (rad * 180 / pi);
}

double distanceEarth(Airport* air1, Airport* air2) {
	return distanceEarth(air1->getLong(), air1->getLat(), air2->getLong(), air2->getLat());
}

	/**
	 * Returns the distance between two points on the Earth.
	 * Direct translation from http://en.wikipedia.org/wiki/Haversine_formula
	 * @file main.cpp The best file
	 * @param lat1d Latitude of the first point in degrees
	 * @param lon1d Longitude of the first point in degrees
	 * @param lat2d Latitude of the second point in degrees
	 * @param lon2d Longitude of the second point in degrees
	 * @return The distance between the two points in miles
	 */
	double distanceEarth(double lat1d, double lon1d, double lat2d, double lon2d) {
		double lat1r, lon1r, lat2r, lon2r, u, v;
		lat1r = deg2rad(lat1d);
		lon1r = deg2rad(lon1d);
		lat2r = deg2rad(lat2d);
		lon2r = deg2rad(lon2d);
		u = sin((lat2r - lat1r)/2);
		v = sin((lon2r - lon1r)/2);
		return 2.0/1.609 * earthRadiusKm * asin(sqrt(u * u + cos(lat1r) * cos(lat2r) * v * v));
	}

void simpleSortTotal(Slist* list, Airport* reference, int low, int high) {
	if(low < high) {
		int m = low+(high-low)/2;
		simpleSortTotal(list, reference, low, m);
		simpleSortTotal(list, reference, m+1, high);
		merge(list, reference, low, m, high);
	}
	//cout << list->rigorSize() << endl;
}

void merge(Slist* list, Airport* reference, int low, int m, int high)
{
	int i, j, k;
	int n1 = m - low + 1;
	int n2 =  high - m;

	Slist* L = list->subList(low, n1);
	//cout << "Sublist: " << low << " " << n1 << endl;
	//cout << L->toString() << endl;
	Slist* R = list->subList(m+1, n2);

	i = 0;
	j = 0;
	k = low;

	while (i < n1 && j < n2)
	{
		if (((Airport*)(L->get(i)))->compare((Airport*)(R->get(j)), reference) > 0)
		{
			list->set(k, L->get(i)->copyDataAndPtr());
			i++;
		}
		else
		{
			list->set(k, R->get(j)->copyDataAndPtr());
			j++;
		}
		k++;
	}

	while (i < n1)
	{
		list->set(k, L->get(i)->copyDataAndPtr());
		i++;
		k++;
	}

	while (j < n2)
	{
		list->set(k, R->get(j)->copyDataAndPtr());
		j++;
		k++;
	}
}

Slist* distanceList(Slist *list, Airport *ref, double distance) {
	auto * ret = new Slist;
	auto * last = dynamic_cast<Airport *>(list->head->ptr);
	int index = 0;
	int size = list->size();
	while(index < size && last != nullptr) {
		if(distanceEarth(last, ref) <= distance) {
			ret->add(last->copyData());
		}
		last = dynamic_cast<Airport *>(last->ptr);
		++index;
	}
	return ret;
}

Airport *farthest(Slist *list, Airport *ref) {
	auto * last = dynamic_cast<Airport *>(list->head->ptr);
	int index = 0;
	int size = list->size();
	double farthest = 0;
	Airport* a = nullptr;
	while(index < size-1 && last != nullptr) {
		double d = distanceEarth(last, ref);
		if(d > farthest && last != ref) {
			a = dynamic_cast<Airport *>(last->copyData());
			farthest = d;
		}
		last = dynamic_cast<Airport *>(last->ptr);
		++index;
	}
	return a;
}

void heapSort(Slist *list, int n, Airport *reference) {
	for(int i = n/2-1; i>=0; --i) {
		heapify(list, reference, n, i);
	}
	for(int i = n-1; i>=0; --i) {
		list->swap(0, i);

		heapify(list, reference, i, 0);
	}
}

void heapify(Slist *list, Airport *reference, int n, int i) {
	int largest = i;
	int l = 2*i+1;
	int r = 2*i+2;

	if(l<n && dynamic_cast<Airport*>(list->get(l))->compare(dynamic_cast<Airport *>(list->get(largest)), reference) > 0) {
		largest = l;
	}
	if(r < n && dynamic_cast<Airport*>(list->get(r))->compare(dynamic_cast<Airport *>(list->get(largest)), reference) > 0) {
		largest = r;
	}

	if(largest != i) {
		list->swap(i, largest);
		heapify(list, reference, n, largest);
	}
}

void printCode(Slist *list, Airport* reference) {
	Airport* last = dynamic_cast<Airport *>(list->head->ptr);
	while(last->ptr != nullptr) {
		//cout << s << " " << last->toString(ref) << endl;
		cout << distanceEarth(last, reference) << endl;
		last = dynamic_cast<Airport *>(last->ptr);
		//s++;
	}
}