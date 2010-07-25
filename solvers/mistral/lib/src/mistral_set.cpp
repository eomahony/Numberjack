
/*
  Mistral is a constraint satisfaction and optimisation library
  Copyright (C) 2003-2005  Emmanuel Hebrard
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

  The author can be contacted electronically at 
  ehebrard@cse.unsw.edu.au.
*/

#include <mistral_set.h>
#include <assert.h>
#include <iomanip>
#include <math.h>




/**********************************************
 * Set
 **********************************************/


void showUint(uint n, std::ostream & os){
  uint mask=1;
  while(mask){
    if(mask & n) os<<1;
    else os<<0;
    mask = mask << 1;
  }
}


using namespace Mistral;
using namespace std;

List::List()
{
  absidx = NULL;
  size = 0;
  list_ = NULL;
}

List::~List()
{
  delete [] list_;
  delete [] absidx;
}

void List::init(const int lb, 
		const int ub,
		const int l)		
{
  capacity = ub-lb+1;
  list_ = new int[capacity+1];
  list_[capacity] = NOVAL;
  absidx = new unsigned int[capacity];
  index_ = (absidx - lb);

  for(unsigned int i=0; i<capacity; ++i) 
    {
      index_[lb+i] = i;
      list_[i] = (lb+i);
    }
  size = l;
}

void List::init(BitSet& e, const int l)
{
  const int lb = e.min();
  const int ub = e.max();
  capacity = e.size();
  list_ = new int[capacity+1];
  list_[capacity] = NOVAL;
  absidx = new unsigned int[ub-lb+1];
  index_ = (absidx - lb);

  int i=0, j=lb;
  do {
    index_[j] = i;
    list_[i] = j;
    ++i;
  } while( (j = e.next(j)) != NOVAL );
  size = l;
}


void List::print(std::ostream& o) const 
{
  cout << "[";
  for(unsigned int i=0; i<size; ++i) {
    o << " " << list_[i];
    o.flush();
    assert( index_[list_[i]] == i );
  }
  cout << " ]";
}
 



// MistralSet::MistralSet(const int lb, const int ub, const uint p) 
// {
//   init(lb,ub,p);
// }

// void MistralSet::reinit(const int lb, const int ub, const uint p) 
// {
//   table += neg_words;
//   delete [] table;
//   init(lb, ub, p);
// }


// void MistralSet::init(const int sz, const uint p) 
// {
//   pos_words = sz;
//   neg_words = 0;

//   if( sz ) {
//     table = new uint[pos_words];
//     for(uint i=0; i<pos_words; ++i) 
//       table[i]=p;
//   } else table = NULL;
// }

// void MistralSet::init(const int lb, const int ub, const uint p) 
// {
//   //  pos_words = 0;
//   //  neg_words = (int)(floor((double)(lb  )/(double)size_word_bit));
//   //  pos_words = (int)( ceil((double)(ub+1)/(double)size_word_bit));

//   neg_words = (lb >> EXP);
//   pos_words = (ub >> EXP)+1;

//   table = new uint[pos_words-neg_words];
//   for(uint i=0; i<pos_words-neg_words; ++i) 
//     table[i]=p;
//   table[pos_words-neg_words-1] &= 
//     (p >> (size_word_bit-1-(ub & CACHE)));
//   table[0] &= (p << (lb & CACHE));
//   table -= neg_words;
// }


// MistralSet::MistralSet(const MistralSet& s) 
// {
//   clone( s );
// }

// void MistralSet::clone(const MistralSet& s)
// {
//   neg_words = s.neg_words;
//   pos_words = s.pos_words;
//   table = new uint[pos_words-neg_words];
//   memcpy(table, s.table+neg_words,
// 	 size_word_byte*(pos_words-neg_words));
//   table -= neg_words;
// }

// void MistralSet::pointTo(MistralSet& s)
// {
//   neg_words = s.neg_words;
//   pos_words = s.pos_words;
//   table = s.table;
// }

// void MistralSet::pointTo(unsigned int *t)
// {
//   neg_words = 0;
//   pos_words = 1;
//   table = t;
// }

// void MistralSet::copy(const MistralSet& s) 
// {
//   int i = (pos_words > s.pos_words ? s.pos_words : pos_words);
//   int j = (neg_words < s.neg_words ? s.neg_words : neg_words);
//   if( i>j )
//     memcpy(table+j,s.table+j,size_word_byte*(i-j));
//   else clear();
// }

// MistralSet::~MistralSet() 
// {
//   table += neg_words;
//   delete [] table; 
// }

// string MistralSet::toString()const
// {
//   string s;
//   if( !empty() ) {
//     int last = NOVAL, cur=min(), end = max(), aft;
//     bool flag=false;
//     do{
//       aft = next(cur);

//       if(aft != cur+1 || cur != last+1) {
// 	if( flag )
// 	  s += " ";
// 	s += int2string( cur );
// 	flag = true;
//       } else if(flag) {
// 	s += "..";
// 	flag = false;
//       }
//       last = cur;
//       cur = aft;
//     } while( cur != NOVAL && cur != last );
//   }
//   return s;
// }

// void  MistralSet::print(std::ostream & os) const 
// {
//   print( os, "{,}" );
// }

// void  MistralSet::print(std::ostream & os, 
// 			const char *delim) const
// {
//   os << delim[0];
//   if( !empty() ) {
//     int last = NOVAL, cur=min(), end = max(), aft;
//     bool flag=false;
//     do{
//       aft = next(cur);

//       if(aft != cur+1 || cur != last+1) {
// 	if( flag )
// 	  os << delim[1];
// 	os << cur;
// 	flag = true;
//       } else if(flag) {
// 	os << "..";
// 	flag = false;
//       }
//       last = cur;
//       cur = aft;
//     } while( cur != NOVAL && cur != last );
//   }
//   os << delim[2];
// }

// void  MistralSet::printBits(std::ostream & os) const 
// {
//   os << "[";
//   for(int i=neg_words; i<pos_words; ++i)
//     showUint( table[i], os );
//   os << "]";
// }

