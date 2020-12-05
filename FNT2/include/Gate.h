//
//	George O'Neill @ University of York, 2020
//
//	This file holds methods for gates
//
#ifndef GATES_H
#define GATES_H

namespace fnt{	//	create unique working area

class Gate{	//	object for channels

public:
	Gate( char c, unsigned long long min, unsigned long long max ){	//	default constructor
		type = c;	//	set gate type
		low = min;	//	set gate minimum
		high = max;	//	set gate maximum
	}	//	end default constructor
	~Gate(){}	//	destructor

	unsigned long long min(){
		return low;
	}	//	get gate minimum
	unsigned long long max(){
		return high;
	}	//	get gate maximum
	double percentage( double n ){
		return (double)( high - low ) / n;
	}	//	get percentage covered
	bool passes( double e ){
		return ( (double)low <= e && (double)high >= e );
	}	//	check if number is inside gate for calibrated entries		
	bool passes( unsigned long long e ){
		return ( low <= e && high >= e );
	}	//	check if number is inside gate for raw entries or time

private:
	char type;	//	gate type
	unsigned long long low, high;	//	gate minimum, maximum

};

}	//	end namespace fnt
#endif