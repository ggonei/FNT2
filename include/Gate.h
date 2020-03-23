//
//	George O'Neill @ University of York, 2020/02/13
//
//	This file holds methods for gates
//
#ifndef GATES_H
#define GATES_H
#include <iostream>

namespace fnt {	//	create unique working area

class Gate {	//	object for channels

	public:
		Gate(char c, ULong64_t min, ULong64_t max) {	//	default constructor
			type = c;	//	set gate type
			low = min;	//	set gate minimum
			high = max;	//	set gate maximum
		}	//	end default constructor
		~Gate() {}	//	destructor
		
		ULong64_t min() { return low; }	//	get gate minimum
		ULong64_t max() { return high; }	//	get gate maximum
		bool passes( ULong64_t e ) { return (min() <= e && max() >= e); }	//	check if energy is inside gate( energy )		

	private:
		char type;	//	gate type
		ULong64_t low, high;	//	gate minimum, maximum

	};

}	//	end namespace fnt
#endif