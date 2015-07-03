
/*
  Mistral 2.0 is a constraint satisfaction and optimisation library
  Copyright (C) 2009  Emmanuel Hebrard
  
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

  The author can be contacted electronically at emmanuel.hebrard@gmail.com.
*/

#include <math.h>

#include <mistral_structure.hpp>
#include <mistral_variable.hpp>
#include <mistral_global.hpp>



void Mistral::IntStack::initialise(IntStack& shared, const int sz) {
  index_capacity = shared.index_capacity;
  
  list_capacity = sz;
  list_ = new int[list_capacity];
  
  start_ = shared.start_;
  index_ = shared.index_;
  
  size = 0;
}

void Mistral::IntStack::initialise(const int lb, const int ub, const int sz, const bool full)
{
  
  index_capacity = ub-lb+1;
  
  list_capacity = sz;
  list_ = new int[list_capacity];
  start_ = new unsigned int[index_capacity];
  index_ = start_ - lb;      
  
  for(int i=lb; i<=ub; ++i) 
    {
      index_[i] = i-lb;
      if(i-lb < (int)list_capacity) list_[i-lb] = i;
    }
  
  size = (full ? sz : 0);
}

std::string Mistral::IntStack::to_str() const {
    int min =  INFTY;
    int max = -INFTY;
  
    if(size) {
      for(unsigned int i=0; i<size; ++i) {
        if(list_[i] < min) min = list_[i];
        if(list_[i] > max) max = list_[i];
      } 
    } else {
      min = max = 0;
    }
  
    BitSet elts(min, max, BitSet::empt);
    for(unsigned int i=0; i<size; ++i) {
      elts.add(list_[i]);
    }
  
    return elts.to_str();
}

std::ostream& Mistral::IntStack::display(std::ostream& os) const {
  int min =  INFTY;
  int max = -INFTY;
  
  if(size) {
    for(unsigned int i=0; i<size; ++i) {
      if(list_[i] < min) min = list_[i];
      if(list_[i] > max) max = list_[i];
    } 
  } else {
    min = max = 0;
  }
  
  BitSet elts(min, max, BitSet::empt);
  for(unsigned int i=0; i<size; ++i) {
    elts.add(list_[i]);
  }
  
  os << elts;

  if(elts.size() != size) {
    std::cout << "ERROR " << elts.size() << " / " << size << std::endl;
    exit(1);
  }

  // os << "(";
  // bool not_empty = (size>0);
  
  // if(not_empty) os << list_[0];
  // for(unsigned int i=1; i<size; ++i)
  // 	os << " " << list_[i];
  // os << ")";
  return os;
}


void Mistral::IntStack::extend_list()
{
  unsigned int increment = ((list_capacity+1) << 1);
  list_capacity += increment;
  
  int* new_list = new int[list_capacity];
  memcpy(new_list, list_, (list_capacity-increment)*sizeof(int));
  
  // for(unsigned int i=0; i<list_capacity-increment; ++i)
  //  	new_list[i] = list_[i];
  
  delete [] list_;
  list_ = new_list;
}

void Mistral::IntStack::extend(const int new_elt)
    {
      int lb = (int)(start_-index_), new_lb = lb;
      int ub = index_capacity+lb-1, new_ub = ub;
      if(new_elt < lb) {
	new_lb = new_elt;
      } else if(new_elt > ub) {
	new_ub = new_elt;
      } else {
	return;
      }
      
      unsigned int new_index_capacity = new_ub-new_lb+1;
      if(new_index_capacity < index_capacity*2) new_index_capacity = index_capacity*2;
      if(new_lb < lb) {
	new_lb = ub-new_index_capacity+1;
      } else {
	new_ub = lb+new_index_capacity-1;
      }



      unsigned int *aux_start = start_;
      start_ = new unsigned int[new_index_capacity];
      std::fill(start_+index_capacity, start_+new_index_capacity, INFTY);
      memcpy(start_+(lb-new_lb), aux_start, index_capacity*sizeof(unsigned int));
      delete [] aux_start;

      index_ = start_ - new_lb;
      int k = 0;
      for(int i=new_lb; i<lb; ++i) {
	index_[i] = size+k++;
	//list_[index_capacity+k++] = i;
      }
      for(int i=ub+1; i<=new_ub; ++i) {
	index_[i] = size+k++;
	//list_[index_capacity+k++] = i;
      }

      index_capacity = new_index_capacity;
    }
    //@}    

    /*!@name Accessors*/
    //@{ 
 int Mistral::IntStack::get_min() const 
    {
      int the_min = INFTY;
      int offset = (index_ - start_);
      unsigned int explored;
      int stop_crit = size+offset;

      for(explored = 0; explored < size 
	    //&& size - explored <= (the_min - offset); 
	    && the_min + (int)explored >= stop_crit; 
	  ++explored) {
	if(list_[explored] < the_min) the_min = list_[explored];
      }
      if(explored < size) {
	while(offset < the_min && index_[offset]>=size) ++offset;
	the_min = offset;
      }



      // int the_min = INFTY;
      // if(size) {
      // 	the_min = list_[0];
      // 	// ratio size / index_capacity
      // 	if(4*size < index_capacity) {
      // 	  for(unsigned int i=1; i<size; ++i)
      // 	    if(list_[i] < the_min) the_min = list_[i];
      // 	} else {
      // 	  int val=(index_ - start_);
      // 	  while( val<the_min && index_[val] >= size )
      // 	    ++val;
      // 	  the_min = val;
      // 	}
      // }


      return the_min;
    }

 int Mistral::IntStack::get_max() const 
    {
      int the_max = -INFTY;
      if(size) {
	the_max = list_[0];
	// ratio size / index_capacity
	if(4*size < index_capacity) {
	  for(unsigned int i=1; i<size; ++i)
	    if(list_[i] > the_max) the_max = list_[i];
	} else {
	  int val=(index_capacity + index_ - start_);
	  while( val>the_max && index_[val] >= size )
	    --val;
	  the_max = val;
	}
      }
      return the_max;
    }


     bool Mistral::IntStack::safe_contain(const int elt) const 
    {
      
      // std::cout << "does " << this << " contain " << elt << " (absolute index: " 
      // 		<< ((unsigned)(elt-(int)(start_-index_))) << "/" << index_capacity ;
      // if((unsigned)(elt-(int)(start_-index_))<index_capacity) {
      // 	std::cout << "list index: " << index_[elt] ;
      // }
      // std::cout << ") => " << (((unsigned)(elt-(int)(start_-index_))<index_capacity && index_[elt]<size) ? "yep" : "no") << "\n";

      return ((unsigned)(elt-(int)(start_-index_))<index_capacity && index_[elt]<size);
    } 
 
     bool Mistral::IntStack::contain(const int elt) const 
    {
      return index_[elt]<size;
    } 
	
    int Mistral::IntStack::get_index(const int elt) const 
   {
     return index_[elt];
   } 
  
     bool Mistral::IntStack::empty() const 
    {
      return !size;
    } 

     int Mistral::IntStack::next(const int elt) const
    {
      unsigned int idx = index_[elt]+1;
      return (idx < size ? list_[idx] : elt);
    }
     int Mistral::IntStack::prev(const int elt) const
    {
      int idx = index_[elt]-1;
      return (idx >= 0 ? list_[idx] : elt);
    }
    
     int Mistral::IntStack::operator[](const unsigned int idx) const
    {
      return list_[idx];
    }

     int& Mistral::IntStack::operator[](const unsigned int idx)
    {
      return list_[idx];
    }
    //@}

    /*!@name List Manipulation*/
    //@{
     int* Mistral::IntStack::begin() 
    {
      return list_;
    }

     int* Mistral::IntStack::end() 
    {
      return &(list_[size]);
    }

     int* Mistral::IntStack::end_mem() 
    {
      return list_+list_capacity;
    }

     void Mistral::IntStack::fill()
    {
      size = list_capacity;
    }

     void Mistral::IntStack::clear()
    {
      size = 0;
    }
  
     void Mistral::IntStack::set_to(const int elt)
    {
      size=1;
      index_[*list_] = index_[elt];
      list_[index_[elt]] = *list_;
      *list_ = elt;
      index_[elt] = 0;
    }

     void Mistral::IntStack::remove(const int elt)
    {
      --size;
      index_[list_[size]] = index_[elt];
      list_[index_[elt]] = list_[size];
      list_[size] = elt;
      index_[elt] = size;
    }

     int Mistral::IntStack::next()
    {
      return list_[size];
    }

     int Mistral::IntStack::pop()
    {
      return list_[--size];
    }

     int Mistral::IntStack::pop_head()
    {
      --size;
      index_[list_[size]] = 0;
      const int elt = *list_;
      *list_ = list_[size];
      list_[size] = elt;
      index_[elt] = size;
      return elt;
    }

     int Mistral::IntStack::head() const
    {
      return *list_;
    }
    
     int Mistral::IntStack::back() const
    {
      return list_[size-1];
    }

     void Mistral::IntStack::init_add(const int elt)
    {
      if(size == list_capacity) extend_list();
      if((unsigned)(elt-(int)(start_-index_)) >= index_capacity) {
	extend(elt);
      }
      list_[size] = elt;
      index_[elt] = size;
      ++size;
    }

     void Mistral::IntStack::add(const int elt)
    {
      index_[list_[size]] = index_[elt];
      list_[index_[elt]] = list_[size];
      list_[size] = elt;
      index_[elt] = size;
      ++size;
    }
	
    void Mistral::IntStack::set(const int elt, const int idx)
   {
     index_[list_[idx]] = index_[elt];
     list_[index_[elt]] = list_[idx];
     list_[idx] = elt;
     index_[elt] = idx;
   }

     void Mistral::IntStack::safe_add(const int elt)
    {
      if(!safe_contain(elt)) extend(elt);
      add(elt);
    }

    // create a new element that can potentially be outside the bounds
     void Mistral::IntStack::create(const int elt)
    {
      extend(elt);
      add(elt);
    }

     void Mistral::IntStack::ordered_add(const int elt)
    {
      // the first non-element goes where elt was
      index_[list_[size]] = index_[elt];
      list_[index_[elt]] = list_[size];

      int idx = size;
      while(idx && list_[idx-1] > elt) { // push every values greater than elt above elt
	list_[idx] = list_[idx-1];
	index_[list_[idx-1]] = idx;
	--idx;
      }

      list_[idx] = elt;
      index_[elt] = idx;
      ++size;
    }


     void Mistral::IntStack::revert_to(const int level)
    {
      size = level;
    }

     void Mistral::IntStack::index()
    {
      for(unsigned int i=0; i<list_capacity; ++i)
	index_[list_[i]] = i;
    }
    //@}




std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::Explanation& x) {
  return x.display(os);
}

std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::Interval& x) {
  return x.display(os);
}

std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::BiInterval& x) {
  return x.display(os);
}

std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::IntStack& x) {
  return x.display(os);
}

std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::Graph& x) {
  return x.display(os);
}

std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::Queue& x) {
  return x.display(os);
}

std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::MultiSet& x) {
  return x.display(os);
}

// std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::VariableQueue& x) {
//   return x.display(os);
// }

std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::Explanation* x) {
  return (x ? x->display(os) : os);
}

std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::IntStack* x) {
  return (x ? x->display(os) : os);
}

std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::Graph* x) {
  return (x ? x->display(os) : os);
}

std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::Queue* x) {
  return (x ? x->display(os) : os);
}

std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::MultiSet* x) {
  return (x ? x->display(os) : os);
}

// std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::VariableQueue* x) {
//   return x->display(os);
// }




int Mistral::__modulo_fct__(const int x, const int m) {
  int mod = x%m;
  if(mod && (mod<0) != (m<0))  mod += m;
  return mod;
}


// /*
//   [HERE WE ASSUME 0<k and 0<=a<=b]

// ================================
// |  min([c, d]) such that       |
// |  ([a, b] % k) = [c, d]       |
// |                              |
// |  if k>b: a                   |
// |                              |
// |  if a<=k<=b: 0               |
// |                              |
// |  so we can assume that k < a |
// |                              |
// |  if (a%k + b - a) >= k: 0    |
// |  else: a%k                   |
// ================================

// [IF k>0 THEN:]

// ================================
// |  min([c, d]) such that       |
// |  ([a, b] % k) = [c, d]       |
// |                              |
// |  if k>b: a                   |
// |                              |
// |  if a<=k<=b: 0               |
// |                              |
// |  so we can assume that k < a |
// |                              |
// |  if (a%k + b - a) >= k: 0    |
// |  else: a%k                   |
// ================================


// */
// int min_modulo(const int a, const int b, const int k) {
//   int value = a;
//   if(k<=b) {
//     if(k>=a) value = 0;
//     else {
//       int mod = __modulo_fct__(a,k);
//       if((mod + b - a) >= k) value = 0;
//       else value = mod;
//     }
//   }
//   return value;
// }
// /*
// ================================
// |  max([c, d]) such that       |
// |  ([a, b] % k) = [c, d]       |
// |                              |
// |  if k>b: b                   |
// |                              |
// |  if a<=k<=b: k-1             |
// |                              |
// |  so we can assume that k < a |
// |                              |
// |  if (b%k - b + a) < 0: k-1   |
// |  else: b%k                   |
// ================================
// */
// int max_modulo(const int a, const int b, const int k) {
//   int value = b;
//   if(k<=b) {
//     if(k>=a) value = k-1;
//     else {
//       int mod = __modulo_fct__(b,k);
//       if((mod - b + a) < 0) value = k-1;
//       else value = mod;
//     }
//   }
//   return value;
// }
// /*
// ================================
// |  min([a, b]) such that       |
// |  ([a, b] % k) = [c, d]       |
// |                              |
// |  if(a%k < c)                 |
// |    a + (c - a%k)             |
// |  if(a%k > d)                 |
// |    a + (d - a%k) + k         |
// |  otherwise: a                |
// ================================
// */
// int min_antimodulo(const int a, const int c, const int d, const int k) {
//   int value = a, mod = __modulo_fct__(a,k);
//   if(mod < c) 
//     value = a + c - mod;
//   else if(mod > d) 
//     value = a + d - mod + k;
//   return value;
// }
// /*
// ================================
// |  max([a, b]) such that       |
// |  ([a, b] % k) = [c, d]       |
// |                              |
// |  if(b%k < c)                 |
// |    b - (b%k - c) - k         |
// |  if(b%k > d)                 |
// |    b - (b%k - d)             |
// |  otherwise: b                |
// ================================
// */
// int max_antimodulo(const int b, const int c, const int d, const int k) {
//   int value = b, mod = __modulo_fct__(b,k);
//   if(mod < c) 
//     value = b - mod + c - k;
//   else if(mod > d) 
//     value = b - mod + d;
//   return value;
// }

Mistral::Interval::Interval() {min=+INFTY; max=-INFTY;}
Mistral::Interval::Interval(const BiInterval b) {
  if(b.positive.empty()) max = (b.zero ? 0 : b.negative.max);
  else max = b.positive.max;

  if(b.negative.empty()) min = (b.zero ? 0 : b.positive.min);
  else min = b.negative.min;
}
Mistral::Interval::Interval(const int _min, const int _max) {min=_min; max=_max;}
Mistral::Interval::~Interval() {}

Mistral::Interval Mistral::Interval::get_union(Mistral::Interval dom) {

  //std::cout << "[" << min << "," << max << "] U [" << dom.min << "," << dom.max << "] = [";

  Interval I((min<=dom.min ? min : dom.min), (max>=dom.max ? max : dom.max));

  //std::cout << I.min << "," << I.max << "]\n";

  return I;
}
  
bool Mistral::Interval::contain(const int x) const {
  return (min <= x && x <= max);
}

bool Mistral::Interval::empty() const { return min>max; }

Mistral::Interval Mistral::Interval::operator-() const {
  Interval I(-max, -min);
  return I;
}

// the interval I such that this%mod = I
Mistral::Interval Mistral::Interval::positive_modulo(const int mod) {
  // the value of the modulo increase with higher values, unless 
  // unless they are distant enough to go back to the next cycle
  Interval I;

  I.min = min%mod;
  I.max = max%mod;

  if(min + mod - I.min <= max) I.min = 0; // can we reach the next 0?
  if(max - I.max - 1 >= min) I.max = mod-1; // can we reach the next mod-1?

  return I;
}

// the interval I such that this%mod = I
Mistral::Interval Mistral::Interval::target_positive_modulo(const int mod, const Interval J) {
  // the value of the modulo increase with higher values, unless 
  // unless they are distant enough to go back to the next cycle
  Interval I;

  //std::cout << "compute pos modulo " << (*this) << " % " << mod << " inter " << J << std::endl; 

  int I_min = I.min = min%mod;
  int I_max = I.max = max%mod;

  if(min + mod - I_min <= max) I.min = 0; // can we reach the next 0?
  if(max - I_max - 1 >= min) I.max = mod-1; // can we reach the next mod-1?


  //std::cout << "[" << I.min << ".." << I_max << "] u [" << I_min << ".." << I.max << "]" << std::endl;


  bool hollow = (max-min+1 < mod && I.min == 0 && I.max == mod-1);

  // do we need to update I.min?
  if(J.min > I.min) { 
    if(
       // is it a hollow interval [0..I_max] u [I_min..mod-1]?
       //max-min+1 < mod && I.min == 0 && I.max == mod-1 && 
       hollow &&
       // is min(J) in the hole?
       J.min > I_max && J.min < I_min
       ) {
      I.min = I_min; // then min(I) takes the next available value
    } else {
      I.min = J.min;
    }
  }

  // std::cout << " max-min+1 < mod = "  
  // 	    << (max-min+1 < mod)
  // 	    << "\n I.min == 0 = " 
  // 	    << (I.min == 0) 
  // 	    << "\n I.max == mod-1 = " 
  // 	    << (I.max == mod-1) 
  // 	    << "\n J.max > I_max = "  
  // 	    << (J.max > I_max)  
  // 	    << "\n J.max < I_min = "
  // 	    << (J.max < I_min) << std::endl;

  // do we need to update I.max?
  if(J.max < I.max) { 
    if(
       // is it a hollow interval [0..I_max] u [I_min..mod-1]?
       //max-min+1 < mod && I.min == 0 && I.max == mod-1 && 
       hollow &&
       // is max(J) in the hole?
       J.max > I_max && J.max < I_min
       ) {
      I.max = I_max; // then max(I) takes the next available value
    } else {
      I.max = J.max;
    }
  }

  // std::cout << " ---> " << I << std::endl;

  return I;
}


// the interval I such that this%mod = I
Mistral::Interval Mistral::Interval::operator%(const int mod) {
  int modulo = abs(mod);
  Interval J;
  if(min>=0) J = positive_modulo(modulo);
  else if(max < 0) {
    Interval I(-max, -min);
    J = -I.positive_modulo(modulo);
  } else {
    Interval Ipos(0, max);
    Interval Ineg(0, -min);

    Interval pos_mod =  Ipos.positive_modulo(modulo);
    Interval neg_mod = -Ineg.positive_modulo(modulo);
    
    J.min = neg_mod.min;
    J.max = pos_mod.max;
  }

  return J;
}

// the interval I such that this%mod = I
Mistral::Interval Mistral::Interval::target_c_modulo(const int mod, const Interval J) {

  //std::cout << "compute c modulo " << (*this) << " % " << mod << " inter " << J << std::endl; 

  int modulo = abs(mod);
  Interval K;
  if(min>=0) {
    K = target_positive_modulo(modulo, J);
  } else if(max < 0) {
    Interval I(-max, -min);
    K = -I.target_positive_modulo(modulo, -J);
  } else {
    Interval Ipos(0, max);
    Interval Ineg(0, -min);

    Interval pos_mod =  Ipos.target_positive_modulo(modulo,  J);
    Interval neg_mod = -Ineg.target_positive_modulo(modulo, -J);
    
    K.min = neg_mod.min;
    K.max = pos_mod.max;
  }

  return K;
}

// the interval I such that this%mod = I
Mistral::Interval Mistral::Interval::operator_modulo(const int mod) {
  // the value of the modulo increase with higher values, unless 
  // unless they are distant enough to go back to the next cycle
  
  int I_min = __modulo_fct__(min,mod);
  int I_max = __modulo_fct__(max,mod);

  if(mod>0) { // positive modulo, hence I \in [0, INFTY]
    if(min + mod - I_min <= max) I_min = 0; // can we reach the next 0?
    if(max - I_max - 1 >= min) I_max = mod-1; // can we reach the next mod-1?
  } else { // negtive modulo, hence I \in [-INFTY, 0]
    if(min - I_min + 1 <= max) I_min = 1+mod;
    if(max + mod - I_max >= min) I_max = 0;
  }

  return Interval(I_min, I_max);
}


// the interval I such that this%mod = I and such that I is a subinterval of J
Mistral::Interval Mistral::Interval::target_modulo(const int mod, const Interval J) {
  // the value of the modulo increase with higher values, unless 
  // unless they are distant enough to go back to the next cycle
  
  int I_min = __modulo_fct__(min,mod);
  int I_max = __modulo_fct__(max,mod);

  Interval I;
  I.min = I_min;
  I.max = I_max;

  if(mod>0) { // positive modulo, hence I \in [0, INFTY]
    if(min + mod - I_min <= max) I.min = 0; // can we reach the next 0?
    if(max - I_max - 1 >= min) I.max = mod-1; // can we reach the next mod-1?
 
    bool hollow = (max-min+1 < mod && I.min == 0 && I.max == mod-1);
    //std::cout << *this << "%" << mod << " = " << I << " ^ " << J << " => " ;

    // do we need to update I.min?
    if(J.min > I.min) { 
      if(
	 // is it a hollow interval [0..I_max] u [I_min..mod-1]?
	 //max-min+1 < mod && I.min == 0 && I.max == mod-1 && 
	 hollow &&
	 // is min(J) in the hole?
	 J.min > I_max && J.min < I_min
	 ) {
	I.min = I_min; // then min(I) takes the next available value
      } else {
	I.min = J.min;
      }
    }

    // do we need to update I.max?
    if(J.max < I.max) { 
      if(
	 // is it a hollow interval [0..I_max] u [I_min..mod-1]?
	 //max-min+1 < mod && I.min == 0 && I.max == mod-1 && 
	 hollow &&
	 // is max(J) in the hole?
	 J.max > I_max && J.max < I_min
	 ) {
	I.max = I_max; // then max(I) takes the next available value
      } else {
	I.max = J.max;
      }
    }

    //std::cout << I << std::endl;


  } else { // negtive modulo, hence I \in [-INFTY, 0]
    if(min - I_min + 1 <= max) I.min = 1+mod;
    if(max + mod - I_max >= min) I.max = 0;


    // std::cout << *this << "%" << mod << " = [" << I.min
    // 	      << ".." << I_max << "] u [" << I_min << ".." << I.max       
    // 	      << "] ^ " << J << " => " ;
    // // std::cout << "min-max-1 < mod: " << (min-max-1) << " < " << mod << " ? " << (min-max-1 < mod) << std::endl;
      
    bool hollow = ((min-max-1 > mod) && I.min == 1+mod && I.max == 0);

    // do we need to update I.min?
    if(J.min > I.min) { 
      if(
	 // is it a hollow interval [1+mod..I_max] u [I_min..0]?
	 //(min-max-1 > mod) && I.min == 1+mod && I.max == 0 && 
	 hollow &&
	 // is min(J) in the hole?
	 (J.min > I_max) && J.min < I_min
	 ) {
	I.min = I_min; // then min(I) takes the next available value
      } else {
	I.min = J.min;
      }
    }


    // std::cout << "(J.max < I.max): " << J.max << " < " << I.max << ": " 
    // 	      << (J.max < I.max) << std::endl
	     

    // do we need to update I.min?
    if(J.max < I.max) { 
      if(
	 // is it a hollow interval [1+mod..I_max] u [I_min..0]?
	 //(min-max-1 > mod) && I.min == 1+mod && I.max == 0 && 
	 hollow &&
	 // is min(J) in the hole?
	 (J.max > I_max) && J.max < I_min
	 ) {
	I.max = I_max; // then min(I) takes the next available value
      } else {
	I.max = J.max;
      }
    }
    
    //  std::cout << I << std::endl;
    
  }
  
  return I; //Interval(I_min, I_max);
}

std::ostream& Mistral::Interval::display(std::ostream& os) const {
  os << "[" << min << "," << max << "]" ;
  return os;
}

std::ostream& Mistral::BiInterval::display(std::ostream& os) const {

  if(!negative.empty()) 
    os << negative;
  if(zero)
    os << "[0]";
  if(!positive.empty()) 
    os << positive;

  return os;
}


void Mistral::Interval::operator+=(const int x) {
  min += x;
  max += x;
}

void Mistral::Interval::operator-=(const int x) {
  min -= x;
  max -= x;
}


Mistral::Interval Mistral::Interval::operator*(const Mistral::Interval arg) {
  BiInterval bi(*this);

  // std::cout << "[" << min << "," << max << "] => [" 
  // 	    << bi.negative.min << "," << bi.negative.max << "|"
  // 	    << (bi.zero ? "0|" : ".|")
  // 	    << bi.positive.min << "," << bi.positive.max << "]\n";
    

  BiInterval bj(arg);

  // std::cout << "[" << arg.min << "," << arg.max << "] => [" 
  // 	    << bj.negative.min << "," << bj.negative.max << "|"
  // 	    << (bj.zero ? "0|" : ".|")
  // 	    << bj.positive.min << "," << bj.positive.max << "]\n";



  BiInterval bk = bi*bj;
  Interval I(bk);


  // std::cout << "[" << I.min << "," << I.max << "] <= [" 
  // 	    << bk.negative.min << "," << bk.negative.max << "|"
  // 	    << (bk.zero ? "0|" : ".|")
  // 	    << bk.positive.min << "," << bk.positive.max << "]\n";



  return I;
}

Mistral::Interval Mistral::Interval::anti_mul(const Mistral::Interval arg) {


  BiInterval bi(*this);


  // std::cout << "[" << min << "," << max << "] => [" 
  // 	    << bi.negative.min << "," << bi.negative.max << "|"
  // 	    << (bi.zero ? "0|" : ".|")
  // 	    << bi.positive.min << "," << bi.positive.max << "]\n";
    

  BiInterval bj(arg);


  // std::cout << "[" << arg.min << "," << arg.max << "] => [" 
  // 	    << bj.negative.min << "," << bj.negative.max << "|"
  // 	    << (bj.zero ? "0|" : ".|")
  // 	    << bj.positive.min << "," << bj.positive.max << "]\n";


  BiInterval bk = bi.anti_mul(bj);
  Interval I(bk);

  // std::cout << "[" << I.min << "," << I.max << "] <= [" 
  // 	    << bk.negative.min << "," << bk.negative.max << "|"
  // 	    << (bk.zero ? "0|" : ".|")
  // 	    << bk.positive.min << "," << bk.positive.max << "]\n";


  return I;
}

Mistral::Interval Mistral::Interval::operator/(const Mistral::Interval arg) {
  BiInterval bi(*this);
  BiInterval bj(arg);
  BiInterval bk = bi/bj;
  return Interval(bk);
}

// Mistral::Interval Mistral::Interval::operator%(const int mod) {
//   BiInterval bi(*this);
//   BiInterval bk = bi%mod;
//   return Interval(bk);
// }


Mistral::NegativeHalfDomain::NegativeHalfDomain() : Interval(INFTY, -INFTY) {}  
Mistral::NegativeHalfDomain::NegativeHalfDomain(const int _min, const int _max) : Interval(_min, _max) {}
Mistral::NegativeHalfDomain::~NegativeHalfDomain() {}

// the interval I such that this*arg = I
Mistral::Interval Mistral::NegativeHalfDomain::operator*(const Mistral::PositiveHalfDomain arg) {
  if(empty() || arg.empty()) return Interval();
  return Interval(min*arg.max, max*arg.min);
}
Mistral::Interval Mistral::NegativeHalfDomain::operator*(const Mistral::NegativeHalfDomain arg) {
  if(empty() || arg.empty()) return Interval();
  return Interval(max*arg.max, min*arg.min);
}
Mistral::Interval Mistral::NegativeHalfDomain::operator*(const int arg) {
  Interval I;
  if(!empty()) {
    if(arg < 0)
      I = Interval(max*arg, min*arg);
    else
      I = Interval(min*arg, max*arg); 
  }
  return I;
}

Mistral::Interval Mistral::NegativeHalfDomain::anti_div_X(const Mistral::PositiveHalfDomain arg) {
  if(empty() || arg.empty() || (!arg.min && !arg.max)) return Interval();
  //(let arg be a positive int)
  // what interval, divided by arg will give this

  // the interval must be negative
  // the minimum value that is in [this] (i.e., >= min) when divided by arg is min*arg+arg-1
  // the maximum value that is in [this] (i.e., <= max) when divided by arg is max*arg;
  return Interval(min*arg.max-arg.max+1, max*arg.min);
}

Mistral::Interval Mistral::NegativeHalfDomain::anti_div_X(const Mistral::NegativeHalfDomain arg) {
  if(empty() || arg.empty() || (!arg.min && !arg.max)) return Interval();
  //(let arg be a negative int)
  // what interval, divided by arg will give this
  
  // the interval must be positive
  // the minimum value that is in [this] (i.e., <= max) when divided by arg is max*arg;
  // the maximum value that is in [this] (i.e., <= max) when divided by arg is min*arg+arg+1;
  return Interval(max*arg.max, min*arg.min-arg.min-1);
}

Mistral::Interval Mistral::NegativeHalfDomain::anti_div_X_pos(const int arg) {
  if(empty() || !arg) return Interval();
  //(let arg be a positive int)
  // what interval, divided by arg will give this

  // the interval must be negative
  // the minimum value that is in [this] (i.e., >= min) when divided by arg is min*arg+arg-1
  // the maximum value that is in [this] (i.e., <= max) when divided by arg is max*arg;
  return Interval(min*arg-arg+1, max*arg);
}

Mistral::Interval Mistral::NegativeHalfDomain::anti_div_X_neg(const int arg) {
  if(empty() || !arg) return Interval();
  //(let arg be a negative int)
  // what interval, divided by arg will give this
  
  // the interval must be positive
  // the minimum value that is in [this] (i.e., <= max) when divided by arg is max*arg;
  // the maximum value that is in [this] (i.e., <= max) when divided by arg is min*arg+arg+1;
  return Interval(max*arg, min*arg-arg-1);
}

Mistral::Interval Mistral::NegativeHalfDomain::anti_div_Y(const Mistral::PositiveHalfDomain arg) {
  if(empty() || arg.empty()) return Interval();

  int ub = arg.min/(min-1)-1;
  int lb = arg.max/max;

  return Interval(lb,ub);
  //return Interval((int)(ceil((double)(arg.max)/(double)(max))), (int)(ceil((double)(arg.min)/(double)(min))));
}

Mistral::Interval Mistral::NegativeHalfDomain::anti_div_Y(const Mistral::NegativeHalfDomain arg) {
  if(empty() || arg.empty()) return Interval();

  int lb = arg.max/(min-1)+1;
  int ub = arg.min/max;

  return Interval(lb,ub);
  //return Interval((int)(floor((double)(arg.max)/(double)(min))), (int)(floor((double)(arg.min)/(double)(max))));
}

Mistral::Interval Mistral::NegativeHalfDomain::anti_div_Y_pos(const int arg) {
  if(empty()) return Interval();
  return Interval((int)(floor((double)(arg)/(double)(max))), (int)(ceil((double)(arg)/(double)(min))));
}

Mistral::Interval Mistral::NegativeHalfDomain::anti_div_Y_neg(const int arg) {
  if(empty()) return Interval();
  return Interval((int)(floor((double)(arg)/(double)(min))), (int)(ceil((double)(arg)/(double)(max))));
}


// the interval I such that this*I = arg
Mistral::Interval Mistral::NegativeHalfDomain::anti_mul(const Mistral::PositiveHalfDomain arg) {
  if(empty() || arg.empty()) return Interval();
  return Interval((int)(ceil((double)min/(double)(arg.min))), (int)(floor((double)(max)/(double)(arg.max))));
}
Mistral::Interval Mistral::NegativeHalfDomain::anti_mul(const Mistral::NegativeHalfDomain arg) {
  if(empty() || arg.empty()) return Interval();
  return Interval((int)(ceil((double)max/(double)(arg.min))), (int)(floor((double)(min)/(double)(arg.max))));
}
Mistral::Interval Mistral::NegativeHalfDomain::anti_mul(const int arg) {
  Interval I;
  if(!empty() || !arg) {
    if(arg<0)
      I = Interval((int)(ceil((double)max/(double)(arg))), (int)(floor((double)min/(double)(arg))));
    else
      I = Interval((int)(ceil((double)min/(double)(arg))), (int)(floor((double)max/(double)(arg))));
  }
  return I;
}

// the interval I such that this/arg = I
Mistral::Interval Mistral::NegativeHalfDomain::operator/(const Mistral::PositiveHalfDomain arg) {
  if(empty() || arg.empty()) return Interval();

  //std::cout << min << "/" << arg.min << " , " << max << "/" << arg.max << std::endl;

  return Interval(min/arg.min, max/arg.max);
}
Mistral::Interval Mistral::NegativeHalfDomain::operator/(const Mistral::NegativeHalfDomain arg) {
  if(empty() || arg.empty()) return Interval();
  return Interval(max/arg.min, min/arg.max);
}
Mistral::Interval Mistral::NegativeHalfDomain::operator/(const int arg) {
  Interval I;
  if(!empty() || !arg) {
    if(arg<0)
      I = Interval(max/arg, min/arg);
    else
      I = Interval(min/arg, max/arg);
  }
  return I;
}

Mistral::Interval Mistral::NegativeHalfDomain::divided_by(const Mistral::NegativeHalfDomain Y,
							  const Mistral::Variable target) {

  if(empty() || Y.empty() || target.get_max()<0) return Interval();

  PositiveHalfDomain X(-max, -min);
  PositiveHalfDomain Yb(-(Y.max), -(Y.min));
  return X.divided_by(Yb, target);

//   // an upper bound z \in Z is reachable when dividing X by Y iff
//   // there exists a value x \in X and y \in Y such that int(x/y) = z
//   // to check that, we look at the interval [min(x)/z, max(x)/z], if it contain an integer value, then this is y.
//   // otherwise, we need to look for the next higher y, in order to give the highest possible value smaller than z.
//   // i.e., we consider y = int((max(x)/z)+1), hence the candidate value for z is int(max(x) / y) = int(max(x) / int((max(x)/z)+1))
//   // we set z to the highest value in Z such that z <= int(max(x) / int((max(x)/z)+1)) and repeat the process until
//   // we have found a reachable upper or non exists


// #ifdef _DEBUG_INTERVAL_DIV
//   std::cout << "compute the division of X = " << (*this) << " by Y = " << Y << " targeting Z = " << target.get_domain() << std::endl;
//   std::cout << " -> compute the upper bound of Z\n"; 
// #endif

//   int lbX = min;
//   int ubX = max;
  
//   int lbY = Y.min;
//   int ubY = Y.max;

//   // we work only on the positive side of Z
//   int lbTarget = target.get_min();
//   if(lbTarget < 0) lbTarget = 0;

//   int ubTarget = target.get_max();

//   int lb_z = lbTarget;
//   int ub_z = ubTarget;
  
//   int curY = ubY;
//   Interval reachable(ubX/curY, lbX/curY);

// #ifdef _DEBUG_INTERVAL_DIV
//   std::cout << "   -> the currently valid upper bound of Z is " << ubTarget << "\n"; 
//   std::cout << "   -> the interval " << reachable <<  " can be obtained when divinding X = " << (*this) << " by max(Y) = " << curY << std::endl;
// #endif

//   if(lbX/(curY+1) >= (lbX/curY)-1) {
//     reachable.min = lbX/ubY;
// #ifdef _DEBUG_INTERVAL_DIV
//     std::cout << "  -> however, from now on, the intervals are contiguous, so we go directly to the end: " << reachable << std::endl;
// #endif
//   }
  
//   while(!reachable.contain(ubTarget)) {
//     // we do not reach Z's ub with the interval reachable when dividing by Y's lb
// #ifdef _DEBUG_INTERVAL_DIV
//     std::cout << "    -> however, it does not overlap with ub(Z) = " << ubTarget << std::endl;
// #endif
//     if(reachable.min > ubTarget) {
    

//       int y_lb = (int)ceil((double)(min) / (double)(ubTarget));
//       int y_ub = (int)floor((double)(max) / (double)(ubTarget));
// #ifdef _DEBUG_INTERVAL_DIV
//       std::cout << "    -> it is above target, we compute the interval on y that allow to reach " << ubTarget 
// 		<< ": [" << ((double)(min) / (double)(ubTarget)) << "," << ((double)(max) / (double)(ubTarget)) << "]" << std::endl;
// #endif
//       if(y_lb <= y_ub) {
// #ifdef _DEBUG_INTERVAL_DIV
// 	std::cout << "    -> this interval contains an integer (" << y_lb 
// 		  << "), so " << ubTarget << " is consistent, END\n";
// #endif
// 	break;
//       } else {
// 	curY = y_lb;
// 	reachable = Interval(lbX/curY, ubX/curY);
// #ifdef _DEBUG_INTERVAL_DIV
// 	std::cout << "    -> this interval does not contain an integer, so we jump to the next available: " << y_lb 
// 		  << " and recompute the reachable interval: " << reachable << std::endl;
// #endif
// 	if(ubX/(curY+1) >= (lbX/curY)-1) {
// 	  reachable.min = lbX/ubY;
// #ifdef _DEBUG_INTERVAL_DIV
// 	  std::cout << "    -> however, from now on, the intervals are contiguous, so we go directly to the end: " << reachable << std::endl;
// #endif
// 	}
//       }
//     }

//     if(reachable.max < ubTarget) {
// #ifdef _DEBUG_INTERVAL_DIV
//       std::cout << "    -> it is below target, so we decrease the target: " ;//<< std::endl;
// #endif
//       // we are below the target, so we jump to the next target
//       if(reachable.max >= lb_z) {
//   	ubTarget = target.prev(reachable.max+1);
// #ifdef _DEBUG_INTERVAL_DIV
// 	std::cout << ubTarget << " is the new target\n";
// #endif
//       } else {
// #ifdef _DEBUG_INTERVAL_DIV
// 	std::cout << " no more possibilities, we fail\n";
// #endif
//   	ubTarget = reachable.max;
//   	break;
//       }
//     }
//   }

// #ifdef _DEBUG_INTERVAL_DIV
//   std::cout << "end loop on " << reachable << " => use " << ubTarget << " as new ub for Z" << std::endl;
// #endif

//   //ub_z = ubTarget;


// #ifdef _DEBUG_INTERVAL_DIV
//   std::cout << " -> compute the lower bound of Z\n"; 
// #endif

//   //lbX = min;
//   //ubX = max;
  
//   //lbY = Y.min;
//   //ubY = Y.max;

//   //lbTarget = target.get_min();
//   ub_z = target.get_max();
  
//   curY = ubY;
//   reachable = Interval(lbX/curY, ubX/curY);

// #ifdef _DEBUG_INTERVAL_DIV
//   std::cout << "   -> the currently valid lower bound of Z is " << lbTarget << "\n"; 
//   std::cout << "   -> the interval " << reachable <<  " can be obtained when divinding X = " << (*this) << " by max(Y) = " << curY << std::endl;
// #endif

// //   if(ubX/(curY+1) >= (lbX/curY)-1) {
// //     reachable.min = lbX/ubY;
// // #ifdef _DEBUG_INTERVAL_DIV
// //     std::cout << "  -> however, from now on, the intervals are contiguous, so we go directly to the end: " << reachable << std::endl;
// // #endif
// //   }
  
//   while(!reachable.contain(lbTarget)) {
//     // we do not reach Z's lb with the interval reachable when dividing by Y's ub
// #ifdef _DEBUG_INTERVAL_DIV
//     std::cout << "    -> however, it does not overlap with lb(Z) = " << lbTarget << std::endl;
// #endif
//     if(reachable.min < lbTarget) {
    
//       int y_lb = (int)ceil((double)(min) / (double)(lbTarget));
//       int y_ub = (int)floor((double)(max) / (double)(lbTarget));
// #ifdef _DEBUG_INTERVAL_DIV
//       std::cout << "    -> it is below target, we compute the interval on y that allow to reach " << lbTarget 
// 		<< ": [" << ((double)(min) / (double)(lbTarget)) << "," << ((double)(max) / (double)(lbTarget)) << "]" << std::endl;
// #endif
//       if(y_lb <= y_ub) {
// #ifdef _DEBUG_INTERVAL_DIV
// 	std::cout << "    -> this interval contains an integer (" << y_ub 
// 		  << "), so " << lbTarget << " is consistent, END\n";
// #endif
// 	break;
//       } else {
// 	curY = y_ub;
// 	reachable = Interval(lbX/curY, ubX/curY);
// #ifdef _DEBUG_INTERVAL_DIV
// 	std::cout << "    -> this interval does not contain an integer, so we jump to the next available: " << y_ub 
// 		  << " and recompute the reachable interval: " << reachable << std::endl;
// #endif
// // 	if(ubX/(curY+1) >= (lbX/curY)-1) {
// // 	  reachable.min = lbX/ubY;
// // #ifdef _DEBUG_INTERVAL_DIV
// // 	  std::cout << "    -> however, from now on, the intervals are contiguous, so we go directly to the end: " << reachable << std::endl;
// // #endif
// // 	}
//       }
//     }

//     if(reachable.max > lbTarget) {
// #ifdef _DEBUG_INTERVAL_DIV
//       std::cout << "    -> it is above target, so we increase the target: " ;//<< std::endl;
// #endif
//       // we are below the target, so we jump to the next target
//       if(reachable.min <= ub_z) {
//   	lbTarget = target.next(reachable.min-1);
// #ifdef _DEBUG_INTERVAL_DIV
// 	std::cout << lbTarget << " is the new target\n";
// #endif
//       } else {
// #ifdef _DEBUG_INTERVAL_DIV
// 	std::cout << " no more possibilities, we fail\n";
// #endif
//   	lbTarget = reachable.min;
//   	break;
//       }
//     }
//   }

// #ifdef _DEBUG_INTERVAL_DIV
//   std::cout << "end loop on " << reachable << " => use " << lbTarget << " as new lb for Z" << std::endl;
// #endif


//   Interval I(lbTarget, ubTarget);

// #ifdef _DEBUG_INTERVAL_DIV
//   std::cout << " => return " << I << std::endl;
// #endif

//   //exit(1);


//   return I;
}


int get_next_higher_abs(const Mistral::Variable X, const int val, bool pos=true) {
  if(pos) {
    return X.next(val);
  } 
  return -(X.prev(-val));
}

int get_next_lower_abs(const Mistral::Variable X, const int val, bool pos=true) {
  if(pos) {
    return X.prev(val);
  } 
  return -(X.next(-val));
}


Mistral::Interval Mistral::NegativeHalfDomain::divided_by(const Mistral::PositiveHalfDomain Y,
							  const Mistral::Variable target) {
  if(empty() || Y.empty() || target.get_min()>0 ) return Interval();


  PositiveHalfDomain X(-max, -min);
  //PositiveHalfDomain Yb(-(Y.max), -(Y.min));
  //PositiveHalfDomain X = -(*this);
  return -(X.divided_by(Y, target, false));


//   // an upper bound z \in Z is reachable when dividing X by Y iff
//   // there exists a value x \in X and y \in Y such that int(x/y) = z
//   // to check that, we look at the interval [min(x)/z, max(x)/z], if it contain an integer value, then this is y.
//   // otherwise, we need to look for the next higher y, in order to give the highest possible value smaller than z.
//   // i.e., we consider y = int((max(x)/z)+1), hence the candidate value for z is int(max(x) / y) = int(max(x) / int((max(x)/z)+1))
//   // we set z to the highest value in Z such that z <= int(max(x) / int((max(x)/z)+1)) and repeat the process until
//   // we have found a reachable upper or non exists


// #ifdef _DEBUG_INTERVAL_DIV
//   std::cout << "compute the division of X = " << (*this) << " by Y = " << Y << " targeting Z = " << target.get_domain() << std::endl;
//   std::cout << " -> compute the lower bound of Z\n"; 
// #endif

//   int lbX = min;
//   int ubX = max;
  
//   int lbY = Y.min;
//   int ubY = Y.max;

//   int lbTarget = target.get_min();

//   // we work only on the positive side of Z
//   int ubTarget = target.get_max();
//   if(ubTarget > 0) ubTarget = 0;

//   int lb_z = lbTarget;
//   int ub_z = ubTarget;
  
//   int curY = ubY;
//   Interval reachable(ubX/curY, lbX/curY);

// #ifdef _DEBUG_INTERVAL_DIV
//   std::cout << "   -> the currently valid lower bound of Z is " << lbTarget << "\n"; 
//   std::cout << "   -> the interval " << reachable <<  " can be obtained when divinding X = " << (*this) << " by max(Y) = " << curY << std::endl;
// #endif

//   if(ubX/(curY-1) <= (lbX/curY)+1) {
//     reachable.max = lbX/lbY;
// #ifdef _DEBUG_INTERVAL_DIV
//     std::cout << "  -> however, from now on, the intervals are contiguous, so we go directly to the end: " << reachable << std::endl;
// #endif
//   }
  
//   while(!reachable.contain(lbTarget)) {
//     // we do not reach Z's lb with the interval reachable when dividing by Y's lb
// #ifdef _DEBUG_INTERVAL_DIV
//     std::cout << "    -> however, it does not overlap with ub(Z) = " << ubTarget << std::endl;
// #endif
//     if(reachable.max < lbTarget) {
    

//       int y_lb = (int)ceil((double)(max) / (double)(lbTarget));
//       int y_ub = (int)floor((double)(min) / (double)(lbTarget));
// #ifdef _DEBUG_INTERVAL_DIV
//       std::cout << "    -> it is below target, we compute the interval on y that allow to reach " << lbTarget 
// 		<< ": [" << ((double)(max) / (double)(lbTarget)) << "," << ((double)(min) / (double)(lbTarget)) << "]" << std::endl;
// #endif
//       if(y_lb <= y_ub) {
// #ifdef _DEBUG_INTERVAL_DIV
// 	std::cout << "    -> this interval contains an integer (" << y_ub 
// 		  << "), so " << lbTarget << " is consistent, END\n";
// #endif
// 	break;
//       } else {
// 	curY = y_ub;
// 	reachable = Interval(ubX/curY, lbX/curY);
// #ifdef _DEBUG_INTERVAL_DIV
// 	std::cout << "    -> this interval does not contain an integer, so we jump to the next available: " << y_ub 
// 		  << " and recompute the reachable interval: " << reachable << std::endl;
// #endif
// 	if(ubX/(curY-1) <= (lbX/curY)+1) {
// 	  reachable.max = lbX/lbY;
// #ifdef _DEBUG_INTERVAL_DIV
// 	  std::cout << "    -> however, from now on, the intervals are contiguous, so we go directly to the end: " << reachable << std::endl;
// #endif
// 	}
//       }
//     }

//     if(reachable.min > lbTarget) {
// #ifdef _DEBUG_INTERVAL_DIV
//       std::cout << "    -> it is above target, so we increase the target: " ;//<< std::endl;
// #endif
//       // we are below the target, so we jump to the next target
//       if(reachable.min <= ub_z) {
//   	lbTarget = target.next(reachable.min-1);
// #ifdef _DEBUG_INTERVAL_DIV
// 	std::cout << lbTarget << " is the new target\n";
// #endif
//       } else {
// #ifdef _DEBUG_INTERVAL_DIV
// 	std::cout << " no more possibilities, we fail\n";
// #endif
//   	lbTarget = reachable.min;
//   	break;
//       }
//     }
//   }

// #ifdef _DEBUG_INTERVAL_DIV
//   std::cout << "end loop on " << reachable << " => use " << lbTarget << " as new lb for Z" << std::endl;
// #endif

//   //ub_z = ubTarget;


// #ifdef _DEBUG_INTERVAL_DIV
//   std::cout << " -> compute the upper bound of Z\n"; 
// #endif

//   //lbX = min;
//   //ubX = max;
  
//   //lbY = Y.min;
//   //ubY = Y.max;

//   //lbTarget = target.get_min();
//   lb_z = target.get_min();
  
//   curY = lbY;
//   reachable = Interval(ubX/curY, lbX/curY);

// #ifdef _DEBUG_INTERVAL_DIV
//   std::cout << "   -> the currently valid upper bound of Z is " << ubTarget << "\n"; 
//   std::cout << "   -> the interval " << reachable <<  " can be obtained when divinding X = " << (*this) << " by min(Y) = " << curY << std::endl;
// #endif

// //   if(ubX/(curY+1) >= (lbX/curY)-1) {
// //     reachable.min = lbX/ubY;
// // #ifdef _DEBUG_INTERVAL_DIV
// //     std::cout << "  -> however, from now on, the intervals are contiguous, so we go directly to the end: " << reachable << std::endl;
// // #endif
// //   }
  
//   while(!reachable.contain(ubTarget)) {
//     // we do not reach Z's ub with the interval reachable when dividing by Y's lb
// #ifdef _DEBUG_INTERVAL_DIV
//     std::cout << "    -> however, it does not overlap with ub(Z) = " << ubTarget << std::endl;
// #endif
//     if(reachable.min > ubTarget) {
    
//       int y_lb = (int)ceil((double)(max) / (double)(ubTarget));
//       int y_ub = (int)floor((double)(min) / (double)(ubTarget));
// #ifdef _DEBUG_INTERVAL_DIV
//       std::cout << "    -> it is above target, we compute the interval on y that allow to reach " << ubTarget 
// 		<< ": [" << ((double)(max) / (double)(ubTarget)) << "," << ((double)(min) / (double)(ubTarget)) << "]" << std::endl;
// #endif
//       if(y_lb <= y_ub) {
// #ifdef _DEBUG_INTERVAL_DIV
// 	std::cout << "    -> this interval contains an integer (" << y_lb 
// 		  << "), so " << ubTarget << " is consistent, END\n";
// #endif
// 	break;
//       } else {
// 	curY = y_lb;
// 	reachable = Interval(ubX/curY, lbX/curY);
// #ifdef _DEBUG_INTERVAL_DIV
// 	std::cout << "    -> this interval does not contain an integer, so we jump to the next available: " << y_lb 
// 		  << " and recompute the reachable interval: " << reachable << std::endl;
// #endif
// // 	if(ubX/(curY+1) >= (lbX/curY)-1) {
// // 	  reachable.min = lbX/ubY;
// // #ifdef _DEBUG_INTERVAL_DIV
// // 	  std::cout << "    -> however, from now on, the intervals are contiguous, so we go directly to the end: " << reachable << std::endl;
// // #endif
// // 	}
//       }
//     }

//     if(reachable.max < ubTarget) {
// #ifdef _DEBUG_INTERVAL_DIV
//       std::cout << "    -> it is below target, so we decrease the target: " ;//<< std::endl;
// #endif
//       // we are below the target, so we jump to the next target
//       if(reachable.max >= lb_z) {
//   	ubTarget = target.prev(reachable.max+1);
// #ifdef _DEBUG_INTERVAL_DIV
// 	std::cout << ubTarget << " is the new target\n";
// #endif
//       } else {
// #ifdef _DEBUG_INTERVAL_DIV
// 	std::cout << " no more possibilities, we fail\n";
// #endif
//   	ubTarget = reachable.max;
//   	break;
//       }
//     }
//   }

// #ifdef _DEBUG_INTERVAL_DIV
//   std::cout << "end loop on " << reachable << " => use " << ubTarget << " as new ub for Z" << std::endl;
// #endif


//   Interval I(lbTarget, ubTarget);

// #ifdef _DEBUG_INTERVAL_DIV
//   std::cout << " => return " << I << std::endl;
// #endif

//   //exit(1);


//   return I;
}




// the interval I such that I/arg = this
//Mistral::Interval Mistral::NegativeHalfDomain::get_dividand(const Mistral::PositiveHalfDomain arg);
//Mistral::Interval Mistral::NegativeHalfDomain::get_dividand(const Mistral::NegativeHalfDomain arg);
///----> this->operator*(arg) <----///


// the interval I such that arg/I = this
//Mistral::Interval Mistral::NegativeHalfDomain::get_divisor(const Mistral::PositiveHalfDomain arg);
//Mistral::Interval Mistral::NegativeHalfDomain::get_divisor(const Mistral::NegativeHalfDomain arg);
///----> arg->operator/(this) <----///

// // the interval I such that this%mod = I
// Mistral::Interval Mistral::NegativeHalfDomain::operator%(const int mod) {
//   // the value of the modulo increase with higher values, unless 
//   // unless they are distant enough to go back to the next cycle
  
//   int I_min = __modulo_fct__(min,mod);
//   int I_max = __modulo_fct__(max,mod);

//   if(mod>0) { // positive modulo, hence I \in [0, INFTY]
//     if(min + mod - I_min <= max) I_min = 0; // can we reach the next 0?
//     if(max - I_max - 1 >= min) I_max = mod-1; // can we reach the next mod-1?
//   } else { // negtive modulo, hence I \in [-INFTY, 0]
//     if(min - I_min + 1 <= max) I_min = 1+mod;
//     if(max + mod - I_max >= min) I_max = 0
//   }

//   return Interval(I_min, I_max);
// }

// the interval I such that I%mod = this
//Mistral::Interval Mistral::NegativeHalfDomain::anti_modulo(const int mod);


Mistral::PositiveHalfDomain::PositiveHalfDomain() : Interval(INFTY, -INFTY) {}
Mistral::PositiveHalfDomain::PositiveHalfDomain(const int _min, const int _max) : Interval(_min, _max) {}
Mistral::PositiveHalfDomain::~PositiveHalfDomain() {}

Mistral::Interval Mistral::PositiveHalfDomain::operator*(const Mistral::PositiveHalfDomain arg) {
  if(empty() || arg.empty()) return Interval();
  return Interval(min*arg.min, max*arg.max);
}
Mistral::Interval Mistral::PositiveHalfDomain::operator*(const Mistral::NegativeHalfDomain arg) {
  if(empty() || arg.empty()) return Interval();
  return Interval(max*arg.min, min*arg.max);
}
Mistral::Interval Mistral::PositiveHalfDomain::operator*(const int arg) {
  Interval I;
  if(!empty()) {
    if(arg < 0)
      I = Interval(max*arg, min*arg);
    else
      I = Interval(min*arg, max*arg);
  }
  return I;
}


Mistral::Interval Mistral::PositiveHalfDomain::anti_div_X(const Mistral::PositiveHalfDomain arg) {
  if(empty() || arg.empty() || (!arg.min && !arg.max)) return Interval();
  //(let arg be a positive int)
  // what interval, divided by arg will give this
  
  // the interval must be positive
  // the minimum value that is in [this] (i.e., >= min) when divided by arg is min*arg
  // the maximum value that is in [this] (i.e., <= max) when divided by arg is max*arg+arg-1;
  return Interval(min*arg.min, max*arg.max+arg.max-1);
}

Mistral::Interval Mistral::PositiveHalfDomain::anti_div_X(const Mistral::NegativeHalfDomain arg) {
  if(empty() || arg.empty() || (!arg.min && !arg.max)) return Interval();
  //(let arg be a negative int)
  // what interval, divided by arg will give this
  
  // the interval must be negative
  // the minimum value that is in [this] (i.e., <= max) when divided by arg is max*arg+arg+1
  // the maximum value that is in [this] (i.e., <= max) when divided by arg is min*arg;
  return Interval(max*arg.min+arg.min+1, min*arg.max);
}

Mistral::Interval Mistral::PositiveHalfDomain::anti_div_X_pos(const int arg) {
  if(empty() || !arg) return Interval();
  //(let arg be a positive int)
  // what interval, divided by arg will give this
  
  // the interval must be positive
  // the minimum value that is in [this] (i.e., >= min) when divided by arg is min*arg
  // the maximum value that is in [this] (i.e., <= max) when divided by arg is max*arg+arg-1;
  return Interval(min*arg, max*arg+arg-1);
}

Mistral::Interval Mistral::PositiveHalfDomain::anti_div_X_neg(const int arg) {
  if(empty() || !arg) return Interval();
  //(let arg be a negative int)
  // what interval, divided by arg will give this
  
  // the interval must be negative
  // the minimum value that is in [this] (i.e., <= max) when divided by arg is max*arg+arg+1
  // the maximum value that is in [this] (i.e., <= max) when divided by arg is min*arg;
  return Interval(max*arg+arg+1, min*arg);
}


Mistral::Interval Mistral::PositiveHalfDomain::anti_div_Y(const Mistral::PositiveHalfDomain arg) {
  if(empty() || arg.empty()) return Interval();
  
  int lb = arg.min/(max+1)+1;
  int ub = arg.max/min;

  return Interval(lb,ub);

  //return Interval((int)(floor((double)arg.min/(double)(max))), (int)(floor((double)arg.max/(double)(min))));
}

Mistral::Interval Mistral::PositiveHalfDomain::anti_div_Y(const Mistral::NegativeHalfDomain arg) {
  // std::cout << arg << " // ?  = " << (*this) << std::endl;  
  // std::cout << (empty()) << "|" << (arg.empty()) << std::endl;
  if(empty() || arg.empty()) return Interval();

  int ub = arg.max/(max+1)-1;
  int lb = arg.min/min;

  return Interval(lb,ub);
  //  return Interval((int)(ceil((double)arg.min/(double)(min))), (int)(ceil((double)arg.max/(double)(max))));
}

Mistral::Interval Mistral::PositiveHalfDomain::anti_div_Y_pos(const int arg) {
  if(empty()) return Interval();
  return Interval((int)(floor((double)arg/(double)(max))), (int)(ceil((double)arg/(double)(min))));
}

Mistral::Interval Mistral::PositiveHalfDomain::anti_div_Y_neg(const int arg) {
  if(empty()) return Interval();
  return Interval((int)(floor((double)arg/(double)(min))), (int)(ceil((double)arg/(double)(max))));
}


Mistral::Interval Mistral::PositiveHalfDomain::anti_mul(const Mistral::PositiveHalfDomain arg) {
  if(empty() || arg.empty()) return Interval();
  return Interval((int)(ceil((double)min/(double)(arg.max))), (int)(floor((double)max/(double)(arg.min))));
}
Mistral::Interval Mistral::PositiveHalfDomain::anti_mul(const Mistral::NegativeHalfDomain arg) {
  if(empty() || arg.empty()) return Interval();
  return Interval((int)(ceil((double)max/(double)(arg.max))), (int)(floor((double)min/(double)(arg.min))));
}
Mistral::Interval Mistral::PositiveHalfDomain::anti_mul(const int arg) {
  Interval I;
  if(!empty() || !arg) {
    if(arg<0)
      I = Interval((int)(ceil((double)max/(double)(arg))), (int)(floor((double)min/(double)(arg))));
    else
      I = Interval((int)(ceil((double)min/(double)(arg))), (int)(floor((double)max/(double)(arg))));
  }
  return I;
}

Mistral::Interval Mistral::PositiveHalfDomain::operator/(const Mistral::PositiveHalfDomain arg) {
  if(empty() || arg.empty()) return Interval();
  return Interval(min/arg.max, max/arg.min);
}

Mistral::Interval Mistral::PositiveHalfDomain::operator/(const Mistral::NegativeHalfDomain arg) {
  if(empty() || arg.empty()) return Interval();
  return Interval(max/arg.max, min/arg.min);
}

Mistral::Interval Mistral::PositiveHalfDomain::divided_by(const Mistral::PositiveHalfDomain Y,
							  const Mistral::Variable target,
							  const bool pos) {
  if(empty() || Y.empty() 
     || (pos && target.get_max()<0)
     || (!pos && target.get_min()>0))  return Interval();

  // an upper bound z \in Z is reachable when dividing X by Y iff
  // there exists a value x \in X and y \in Y such that int(x/y) = z
  // to check that, we look at the interval [min(x)/z, max(x)/z], if it contain an integer value, then this is y.
  // otherwise, we need to look for the next higher y, in order to give the highest possible value smaller than z.
  // i.e., we consider y = int((max(x)/z)+1), hence the candidate value for z is int(max(x) / y) = int(max(x) / int((max(x)/z)+1))
  // we set z to the highest value in Z such that z <= int(max(x) / int((max(x)/z)+1)) and repeat the process until
  // we have found a reachable upper or non exists


#ifdef _DEBUG_INTERVAL_DIV
  std::cout << "compute the division of X = " << (*this) << " by Y = " << Y << " targeting Z = " << target.get_domain() << std::endl;
  std::cout << " -> compute the upper bound of Z\n"; 
#endif

  int lbX = min;
  int ubX = max;
  
  int lbY = Y.min;
  int ubY = Y.max;

  // we work only on the positive side of Z
  int lbTarget = (pos ? target.get_min() : -(target.get_max()));
  if(lbTarget < 0) lbTarget = 0;

  int ubTarget = (pos ? target.get_max() : -(target.get_min()));

  int lb_z = lbTarget;
  int ub_z = ubTarget;
  
  int curY = lbY;
  Interval reachable(lbX/curY, ubX/curY);

#ifdef _DEBUG_INTERVAL_DIV
  std::cout << "   -> the currently valid upper bound of Z is " << ubTarget << "\n"; 
  std::cout << "   -> the interval " << reachable <<  " can be obtained when divinding X = " << (*this) << " by min(Y) = " << curY << std::endl;
#endif

  if(ubX/(curY+1) >= (lbX/curY)-1) {
    reachable.min = lbX/ubY;
#ifdef _DEBUG_INTERVAL_DIV
    std::cout << "  -> however, from now on, the intervals are contiguous, so we go directly to the end: " << reachable << std::endl;
#endif
  }
  
  while(!reachable.contain(ubTarget)) {
    // we do not reach Z's ub with the interval reachable when dividing by Y's lb
#ifdef _DEBUG_INTERVAL_DIV
    std::cout << "    -> however, it does not overlap with ub(Z) = " << ubTarget << std::endl;
#endif

    if(reachable.min > ubTarget) {

      
      // int y_lb = (int)ceil((double)(min) / (double)(ubTarget));
      // int y_ub = (int)floor((double)(max) / (double)(ubTarget));
      int y_lb = (ubTarget ? min / ubTarget : min+1);
      int y_ub = (ubTarget ? max / ubTarget : INFTY);
#ifdef _DEBUG_INTERVAL_DIV
      std::cout << "    -> it is above target, we compute the interval on y that allow to reach " << ubTarget 
	//<< ": [" << ((double)(min) / (double)(ubTarget)) << "," << ((double)(max) / (double)(ubTarget)) << "]" << std::endl;
		<< ": [" << y_lb << "," << y_ub << "]" << std::endl
		<< "       then maps it back to Z: [ " << (y_lb ? max/y_lb : INFTY) << "," << (y_ub ? min/y_ub : INFTY) << "]" << std::endl; 


#endif
      if((!y_lb || (max/y_lb >= ubTarget)) 
	 && y_ub && ubTarget >= min/y_ub) {
#ifdef _DEBUG_INTERVAL_DIV
	std::cout << "    -> this interval contains " << ubTarget 
		  << " so it is consistent, END\n";
#endif
	break;
      } else {
//       if(y_lb <= y_ub) {
// #ifdef _DEBUG_INTERVAL_DIV
// 	std::cout << "    -> this interval contains an integer (" << y_lb 
// 		  << "), so " << ubTarget << " is consistent, END\n";
// #endif
// 	break;
//       } else {
	curY = y_ub+1;

	if(curY > ubY) {
	  ubTarget = -INFTY;
	  break;
	}

	
#ifdef _DEBUG_INTERVAL_DIV
	std::cout << "    -> this interval does not contain " << ubTarget << ", so we jump to the next available: " << curY << std::endl;
#endif
	reachable = Interval(lbX/curY, ubX/curY);
#ifdef _DEBUG_INTERVAL_DIV
	std::cout << "       and recompute the reachable interval: " << reachable << std::endl;
#endif
	if(ubX/(curY+1) >= (lbX/curY)-1) {
	  reachable.min = lbX/ubY;
#ifdef _DEBUG_INTERVAL_DIV
	  std::cout << "    -> however, from now on, the intervals are contiguous, so we go directly to the end: " << reachable << std::endl;
#endif
	}
      }
    }

    if(reachable.max < ubTarget) {
#ifdef _DEBUG_INTERVAL_DIV
      std::cout << "    -> it is below target, so we decrease the target: " ;//<< std::endl;
#endif
      // we are below the target, so we jump to the next target
      if(reachable.max >= lb_z) {
  	//ubTarget = target.prev(reachable.max+1);
	ubTarget = get_next_lower_abs(target, reachable.max+1, pos);
	if(ubTarget < 0) break;
#ifdef _DEBUG_INTERVAL_DIV
	std::cout << ubTarget << " is the new target\n";
#endif
      } else {
#ifdef _DEBUG_INTERVAL_DIV
	std::cout << " no more possibilities, we fail\n";
#endif
  	ubTarget = reachable.max;
  	break;
      }
    }
  }

#ifdef _DEBUG_INTERVAL_DIV
  std::cout << "end loop on " << reachable << " => use " << ubTarget << " as new ub for Z" << std::endl;
#endif

  //ub_z = ubTarget;


#ifdef _DEBUG_INTERVAL_DIV
  std::cout << " -> compute the lower bound of Z\n"; 
#endif

  //lbX = min;
  //ubX = max;
  
  //lbY = Y.min;
  //ubY = Y.max;

  //lbTarget = target.get_min();
  //ub_z = target.get_max();
  
  curY = ubY;
  reachable = Interval(lbX/curY, ubX/curY);

#ifdef _DEBUG_INTERVAL_DIV
  std::cout << "   -> the currently valid lower bound of Z is " << lbTarget << "\n"; 
  std::cout << "   -> the interval " << reachable <<  " can be obtained when divinding X = " << (*this) << " by max(Y) = " << curY << std::endl;
#endif

//   if(ubX/(curY+1) >= (lbX/curY)-1) {
//     reachable.min = lbX/ubY;
// #ifdef _DEBUG_INTERVAL_DIV
//     std::cout << "  -> however, from now on, the intervals are contiguous, so we go directly to the end: " << reachable << std::endl;
// #endif
//   }
  
  while(!reachable.contain(lbTarget)) {
    // we do not reach Z's lb with the interval reachable when dividing by Y's ub
#ifdef _DEBUG_INTERVAL_DIV
    std::cout << "    -> however, it does not overlap with lb(Z) = " << lbTarget << std::endl;
#endif
    if(reachable.min < lbTarget) {
    
      // int y_lb = (int)ceil((double)(min) / (double)(lbTarget));
      // int y_ub = (int)floor((double)(max) / (double)(lbTarget));
      int y_lb = min / lbTarget;
      int y_ub = max / lbTarget;
#ifdef _DEBUG_INTERVAL_DIV
      std::cout << "    -> it is below target, we compute the interval on y that allow to reach " << lbTarget 
	//<< ": [" << ((double)(min) / (double)(lbTarget)) << "," << ((double)(max) / (double)(lbTarget)) << "]" << std::endl;
		<< ": [" << y_lb << "," << y_ub << "]" << std::endl
		<< "       then maps it back to Z: [ " << (y_lb ? max/y_lb : INFTY) << "," << (y_ub ? min/y_ub : INFTY) << "]" << std::endl; 

#endif
      if((!y_lb || (max/y_lb >= lbTarget)) && 
	 y_ub && lbTarget >= min/y_ub) {
#ifdef _DEBUG_INTERVAL_DIV
	std::cout << "    -> this interval contains " << lbTarget
		  << ", so it is consistent, END\n";
#endif
	break;
      } else {
//       if(y_lb <= y_ub) {
// #ifdef _DEBUG_INTERVAL_DIV
// 	std::cout << "    -> this interval contains an integer (" << y_ub 
// 		  << "), so " << lbTarget << " is consistent, END\n";
// #endif
// 	break;
//       } else {
	curY = y_lb;
#ifdef _DEBUG_INTERVAL_DIV
	std::cout << "    -> this interval does not contain " << lbTarget << ", so we jump to the next available: " << curY << std::endl; 
#endif

	if(curY < lbY) {
	  lbTarget = INFTY;
	  break;
	}
	reachable = Interval(lbX/curY, ubX/curY);

#ifdef _DEBUG_INTERVAL_DIV
	std::cout << "       and recompute the reachable interval: " << reachable << std::endl;
#endif
// 	if(ubX/(curY+1) >= (lbX/curY)-1) {
// 	  reachable.min = lbX/ubY;
// #ifdef _DEBUG_INTERVAL_DIV
// 	  std::cout << "    -> however, from now on, the intervals are contiguous, so we go directly to the end: " << reachable << std::endl;
// #endif
// 	}
      }
    }

    if(reachable.max > lbTarget) {
#ifdef _DEBUG_INTERVAL_DIV
      std::cout << "    -> it is above target, so we increase the target: " ;//<< std::endl;
#endif
      // we are below the target, so we jump to the next target
      if(reachable.min <= ub_z) {
  	//lbTarget = target.next(reachable.min-1);
	lbTarget = get_next_higher_abs(target, reachable.min-1, pos);

#ifdef _DEBUG_INTERVAL_DIV
	std::cout << lbTarget << " is the new target\n";
#endif
      } else {
#ifdef _DEBUG_INTERVAL_DIV
	std::cout << " no more possibilities, we fail\n";
#endif
  	lbTarget = reachable.min;
  	break;
      }
    }
  }

#ifdef _DEBUG_INTERVAL_DIV
  std::cout << "end loop on " << reachable << " => use " << lbTarget << " as new lb for Z" << std::endl;
#endif


  Interval I(lbTarget, ubTarget);

#ifdef _DEBUG_INTERVAL_DIV
  std::cout << " => return " << I << std::endl;
#endif

  //exit(1);


  return I;
}




Mistral::Interval Mistral::PositiveHalfDomain::divided_by(const Mistral::NegativeHalfDomain Y,
							  const Mistral::Variable target) {
  if(empty() || Y.empty() || target.get_min()>0 ) return Interval();

  //PositiveHalfDomain X(-max, -min);
  PositiveHalfDomain Yb(-(Y.max), -(Y.min));
  //PositiveHalfDomain Yb = -Y;
  return -(divided_by(Yb, target, false));


//   // an upper bound z \in Z is reachable when dividing X by Y iff
//   // there exists a value x \in X and y \in Y such that int(x/y) = z
//   // to check that, we look at the interval [min(x)/z, max(x)/z], if it contain an integer value, then this is y.
//   // otherwise, we need to look for the next higher y, in order to give the highest possible value smaller than z.
//   // i.e., we consider y = int((max(x)/z)+1), hence the candidate value for z is int(max(x) / y) = int(max(x) / int((max(x)/z)+1))
//   // we set z to the highest value in Z such that z <= int(max(x) / int((max(x)/z)+1)) and repeat the process until
//   // we have found a reachable upper or non exists


// #ifdef _DEBUG_INTERVAL_DIV
//   std::cout << "compute the division of X = " << (*this) << " by Y = " << Y << " targeting Z = " << target.get_domain() << std::endl;
//   std::cout << " -> compute the lower bound of Z\n"; 
// #endif

//   int lbX = min;
//   int ubX = max;
  
//   int lbY = Y.min;
//   int ubY = Y.max;

//   int lbTarget = target.get_min();

//   // we work only on the positive side of Z
//   int ubTarget = target.get_max();
//   if(ubTarget > 0) ubTarget = 0;

//   int lb_z = lbTarget;
//   int ub_z = ubTarget;
  
//   int curY = ubY;
//   Interval reachable(ubX/curY, lbX/curY);

// #ifdef _DEBUG_INTERVAL_DIV
//   std::cout << "   -> the currently valid lower bound of Z is " << lbTarget << "\n"; 
//   std::cout << "   -> the interval " << reachable <<  " can be obtained when divinding X = " << (*this) << " by max(Y) = " << curY << std::endl;
// #endif

//   if(ubX/(curY-1) <= (lbX/curY)+1) {
//     reachable.max = lbX/lbY;
// #ifdef _DEBUG_INTERVAL_DIV
//     std::cout << "  -> however, from now on, the intervals are contiguous, so we go directly to the end: " << reachable << std::endl;
// #endif
//   }
  
//   while(!reachable.contain(lbTarget)) {
//     // we do not reach Z's lb with the interval reachable when dividing by Y's lb
// #ifdef _DEBUG_INTERVAL_DIV
//     std::cout << "    -> however, it does not overlap with ub(Z) = " << ubTarget << std::endl;
// #endif
//     if(reachable.max < lbTarget) {
    

//       int y_lb = (int)ceil((double)(max) / (double)(lbTarget));
//       int y_ub = (int)floor((double)(min) / (double)(lbTarget));
// #ifdef _DEBUG_INTERVAL_DIV
//       std::cout << "    -> it is below target, we compute the interval on y that allow to reach " << lbTarget 
// 		<< ": [" << ((double)(max) / (double)(lbTarget)) << "," << ((double)(min) / (double)(lbTarget)) << "]" << std::endl;
// #endif
//       if(y_lb <= y_ub) {
// #ifdef _DEBUG_INTERVAL_DIV
// 	std::cout << "    -> this interval contains an integer (" << y_ub 
// 		  << "), so " << lbTarget << " is consistent, END\n";
// #endif
// 	break;
//       } else {
// 	curY = y_ub;
// 	reachable = Interval(ubX/curY, lbX/curY);
// #ifdef _DEBUG_INTERVAL_DIV
// 	std::cout << "    -> this interval does not contain an integer, so we jump to the next available: " << y_ub 
// 		  << " and recompute the reachable interval: " << reachable << std::endl;
// #endif
// 	if(ubX/(curY-1) <= (lbX/curY)+1) {
// 	  reachable.max = lbX/lbY;
// #ifdef _DEBUG_INTERVAL_DIV
// 	  std::cout << "    -> however, from now on, the intervals are contiguous, so we go directly to the end: " << reachable << std::endl;
// #endif
// 	}
//       }
//     }

//     if(reachable.min > lbTarget) {
// #ifdef _DEBUG_INTERVAL_DIV
//       std::cout << "    -> it is above target, so we increase the target: " ;//<< std::endl;
// #endif
//       // we are below the target, so we jump to the next target
//       if(reachable.min <= ub_z) {
//   	lbTarget = target.next(reachable.min-1);
// #ifdef _DEBUG_INTERVAL_DIV
// 	std::cout << lbTarget << " is the new target\n";
// #endif
//       } else {
// #ifdef _DEBUG_INTERVAL_DIV
// 	std::cout << " no more possibilities, we fail\n";
// #endif
//   	lbTarget = reachable.min;
//   	break;
//       }
//     }
//   }

// #ifdef _DEBUG_INTERVAL_DIV
//   std::cout << "end loop on " << reachable << " => use " << lbTarget << " as new lb for Z" << std::endl;
// #endif

//   //ub_z = ubTarget;


// #ifdef _DEBUG_INTERVAL_DIV
//   std::cout << " -> compute the upper bound of Z\n"; 
// #endif

//   //lbX = min;
//   //ubX = max;
  
//   //lbY = Y.min;
//   //ubY = Y.max;

//   //lbTarget = target.get_min();
//   lb_z = target.get_min();
  
//   curY = lbY;
//   reachable = Interval(ubX/curY, lbX/curY);

// #ifdef _DEBUG_INTERVAL_DIV
//   std::cout << "   -> the currently valid upper bound of Z is " << ubTarget << "\n"; 
//   std::cout << "   -> the interval " << reachable <<  " can be obtained when divinding X = " << (*this) << " by min(Y) = " << curY << std::endl;
// #endif

// //   if(ubX/(curY+1) >= (lbX/curY)-1) {
// //     reachable.min = lbX/ubY;
// // #ifdef _DEBUG_INTERVAL_DIV
// //     std::cout << "  -> however, from now on, the intervals are contiguous, so we go directly to the end: " << reachable << std::endl;
// // #endif
// //   }
  
//   while(!reachable.contain(ubTarget)) {
//     // we do not reach Z's ub with the interval reachable when dividing by Y's lb
// #ifdef _DEBUG_INTERVAL_DIV
//     std::cout << "    -> however, it does not overlap with ub(Z) = " << ubTarget << std::endl;
// #endif
//     if(reachable.min > ubTarget) {
    
//       int y_lb = (int)ceil((double)(max) / (double)(ubTarget));
//       int y_ub = (int)floor((double)(min) / (double)(ubTarget));
// #ifdef _DEBUG_INTERVAL_DIV
//       std::cout << "    -> it is above target, we compute the interval on y that allow to reach " << ubTarget 
// 		<< ": [" << ((double)(max) / (double)(ubTarget)) << "," << ((double)(min) / (double)(ubTarget)) << "]" << std::endl;
// #endif
//       if(y_lb <= y_ub) {
// #ifdef _DEBUG_INTERVAL_DIV
// 	std::cout << "    -> this interval contains an integer (" << y_lb 
// 		  << "), so " << ubTarget << " is consistent, END\n";
// #endif
// 	break;
//       } else {
// 	curY = y_lb;
// 	reachable = Interval(ubX/curY, lbX/curY);
// #ifdef _DEBUG_INTERVAL_DIV
// 	std::cout << "    -> this interval does not contain an integer, so we jump to the next available: " << y_lb 
// 		  << " and recompute the reachable interval: " << reachable << std::endl;
// #endif
// // 	if(ubX/(curY+1) >= (lbX/curY)-1) {
// // 	  reachable.min = lbX/ubY;
// // #ifdef _DEBUG_INTERVAL_DIV
// // 	  std::cout << "    -> however, from now on, the intervals are contiguous, so we go directly to the end: " << reachable << std::endl;
// // #endif
// // 	}
//       }
//     }

//     if(reachable.max < ubTarget) {
// #ifdef _DEBUG_INTERVAL_DIV
//       std::cout << "    -> it is below target, so we decrease the target: " ;//<< std::endl;
// #endif
//       // we are below the target, so we jump to the next target
//       if(reachable.max >= lb_z) {
//   	//ubTarget = target.prev(reachable.max+1);
// 	ubTarget = get_next_lower_abs(target, reachable.max+1, pos);
// #ifdef _DEBUG_INTERVAL_DIV
// 	std::cout << ubTarget << " is the new target\n";
// #endif
//       } else {
// #ifdef _DEBUG_INTERVAL_DIV
// 	std::cout << " no more possibilities, we fail\n";
// #endif
//   	ubTarget = reachable.max;
//   	break;
//       }
//     }
//   }

// #ifdef _DEBUG_INTERVAL_DIV
//   std::cout << "end loop on " << reachable << " => use " << ubTarget << " as new ub for Z" << std::endl;
// #endif


//   Interval I(lbTarget, ubTarget);

// #ifdef _DEBUG_INTERVAL_DIV
//   std::cout << " => return " << I << std::endl;
// #endif

//   //exit(1);


//   return I;
}



// Mistral::Interval Mistral::PositiveHalfDomain::operator/(const Mistral::NegativeHalfDomain arg,
// 							 const Mistral::NegativeHalfDomain target) {
//   if(empty() || arg.empty()) return Interval();
//   return Interval(max/arg.max, min/arg.min);
// }
Mistral::Interval Mistral::PositiveHalfDomain::operator/(const int arg) {
  Interval I;
  if(!empty() || !arg) {
    if(arg<0)
      I = Interval(max/arg, min/arg);
    else
      I = Interval(min/arg, max/arg);
  }
  return I;
}

// Mistral::Interval Mistral::PositiveHalfDomain::anti_div_X(const Mistral::PositiveHalfDomain arg);
// Mistral::Interval Mistral::PositiveHalfDomain::anti_div_X(const Mistral::NegativeHalfDomain arg);

// Mistral::Interval Mistral::PositiveHalfDomain::operator%(const int mod) {
//   // the value of the modulo increase with higher values, unless 
//   // unless they are distant enough to go back to the next cycle
  
//   int I_min = __modulo_fct__(min,mod);
//   int I_max = __modulo_fct__(max,mod);

//   if(mod>0) { // positive modulo, hence I \in [0, INFTY]
//     if(min + mod - I_min <= max) I_min = 0; // can we reach the next 0?
//     if(max - I_max - 1 >= min) I_max = mod-1; // can we reach the next mod-1?
//   } else { // negtive modulo, hence I \in [-INFTY, 0]
//     if(min - I_min + 1 <= max) I_min = 1+mod;
//     if(max + mod - I_max >= min) I_max = 0
//   }

//   return Interval(I_min, I_max);
// }

// // // the return value is the interval I such that I%mod = this
// // // 
// // Mistral::Interval Mistral::PositiveHalfDomain::anti_modulo(const int mod) {
// //   int mod_min =  __modulo_fct__(min,mod);

// // }



Mistral::BiInterval::BiInterval(const int n_min, const int n_max, const int p_min, const int p_max, const bool z) {
  positive.min = p_min;
  positive.max = p_max;
  negative.min = n_min;
  negative.max = n_max;
  zero = z;
}
Mistral::BiInterval::BiInterval(const Interval _neg, const Interval _pos, const bool z) {
  positive.min = _pos.min;
  positive.max = _pos.max;
  if(positive.min <  1) positive.min =  1;

  negative.min = _neg.min;
  negative.max = _neg.max;
  if(negative.max > -1) negative.max = -1;

  zero = (z || _pos.contain(0) || _neg.contain(0));
}
Mistral::BiInterval::BiInterval(const Interval _neg, const Interval _pos, const Interval z) {
  positive.min = _pos.min;
  positive.max = _pos.max;
  if(1 > positive.min || z.contain(1)) positive.min =  1;
  if(z.max > positive.max) positive.max = z.max;
  
  negative.min = _neg.min;
  negative.max = _neg.max;
  if(negative.max > -1 || z.contain(-1)) negative.max = -1;
  if(z.min < negative.min) negative.min = z.min;

  zero = (z.contain(0) || _pos.contain(0) || _neg.contain(0));
}
Mistral::BiInterval::BiInterval(const Interval I) {
  initialise(I.min, I.max);  
}

Mistral::BiInterval::BiInterval(const int min, const int max) {
  initialise(min, max);  
}

void Mistral::BiInterval::initialise(const int min, const int max) {

  if(min>=0) {
    positive.min = min+(min==0);
    positive.max = max;

    negative.min = +INFTY;
    negative.max = -INFTY;
    
    zero = (min==0);
  } else if(max<=0) {
    positive.min = +INFTY;
    positive.max = -INFTY;

    negative.min = min;
    negative.max = max-(max==0);
    
    zero = (max==0);
  } else { // min<0 & max>0
    positive.min = 1;
    positive.max = max;

    negative.min = min;
    negative.max = -1;

    zero = true;
  }
}
Mistral::BiInterval::BiInterval() {
  positive.min = +INFTY;
  positive.max = -INFTY;
  
  negative.min = +INFTY;
  negative.max = -INFTY;
  
  zero = false;
}
Mistral::BiInterval::BiInterval(const Variable x) {

  // std::cout << "build a bi-interval from " << x << " in " << x.get_domain() << std::endl;
  // std::cout << x.get_max_neg() << ".." << x.get_min_pos() << std::endl;

  positive.min = x.get_min_pos();
  positive.max = x.get_max();
  
  negative.min = x.get_min();
  negative.max = x.get_max_neg();
  
  zero = x.contain(0);
}
Mistral::BiInterval::~BiInterval() {}


int Mistral::BiInterval::get_min() const {
  int the_min = 0;
  if(negative.empty()) {
    if(!zero) the_min = positive.min;
  } else {
    the_min = negative.min;
  }
  return the_min;
}

int Mistral::BiInterval::get_max() const {
  int the_max = 0;
  if(positive.empty()) {
    if(!zero) the_max = negative.max;
  } else {
    the_max = positive.max;
  }
  return the_max;
}

int Mistral::BiInterval::get_min_abs() const {
  int the_min = 0;
  if(!zero) {
    if(negative.empty()) the_min = positive.min;
    else if(positive.empty()) the_min = -(negative.max);
    else the_min = std::min(positive.min, -(negative.max));
  } 
  return the_min;
}
int Mistral::BiInterval::get_max_abs() const {
  int the_max = 0;
 
  if(negative.empty()) {
    if(!positive.empty()) the_max = positive.max;
  } else if(positive.empty()) the_max = -(negative.min);
  else the_max = std::max(positive.max, -(negative.min));
  
  return the_max;
}

int Mistral::BiInterval::is_hollow() const {
  return (!positive.empty() &&
	  !negative.empty() &&
	  (!zero || negative.max < -1 || positive.min > 1));
}

    
Mistral::BiInterval Mistral::BiInterval::operator*(const int arg) {
  // Interval pos = (arg<0 ? negative : positive) * arg;
  // Interval neg = (arg<0 ? positive : negative) * arg;
  Interval pos, neg;
  if(arg < 0) {
    pos = negative*arg;
    neg = positive*arg;
  } else {
    neg = negative*arg;
    pos = positive*arg;
  }


#ifdef _DEBUG_INTERVALS
  if(_DEBUG_INTERVALS) {
  display(std::cout) ;
  std::cout << " * " << arg << " = ";
  }
#endif

  BiInterval I(neg, pos, zero || !arg);

#ifdef _DEBUG_INTERVALS
  if(_DEBUG_INTERVALS) {
  I.display(std::cout);
  std::cout << std::endl << neg << std::endl << pos << std::endl;
  }
#endif

  return I;
}


Mistral::BiInterval Mistral::BiInterval::anti_div_X(const int arg) {

#ifdef _DEBUG_INTERVALS
  if(_DEBUG_INTERVALS) {
  std::cout << "compute anti div : [min?,max?] / " << arg << " = ";
  display(std::cout);
  std::cout << std::endl;
  }
#endif 

  BiInterval I;

  if(arg) {

    Interval pos, neg, z;
    
    if(arg > 0) {
      pos = positive.anti_div_X_pos(arg);
      neg = negative.anti_div_X_pos(arg);
    } else {
      neg = positive.anti_div_X_neg(arg);
      pos = negative.anti_div_X_neg(arg);
    }

    if(zero) {
      if(arg<0)
	z = Interval(1+arg, -arg-1);
      else
	z = Interval(1-arg,  arg-1);
    }

   
    I = BiInterval(neg, pos, z);

  }

  return I;
}


// X / Y = Z
// X = Z * Y
Mistral::BiInterval Mistral::BiInterval::anti_div_X(const Mistral::BiInterval arg) {

#ifdef _DEBUG_INTERVALS
  if(_DEBUG_INTERVALS) {
  std::cout << "compute anti div : [min?,max?] / " << arg << " = ";
  display(std::cout);
  std::cout << std::endl;
  }
#endif 

  BiInterval I;

  Interval pospos = positive.anti_div_X(arg.positive);

  //std::cout << "pospos: [" << pospos.min << "," << pospos.max << "]\n";

  Interval posneg = positive.anti_div_X(arg.negative);

  //std::cout << "posneg: [" << posneg.min << "," << posneg.max << "]\n";

  Interval negpos = negative.anti_div_X(arg.positive);

  //std::cout << "negpos: [" << negpos.min << "," << negpos.max << "]\n";

  Interval negneg = negative.anti_div_X(arg.negative);

  //std::cout << "negneg: [" << negneg.min << "," << negneg.max << "]\n";

  Interval z;

  if(zero) {
    int amin = -arg.get_min();
    int amax =  arg.get_max();
    int maxabs = std::max(amin, amax);
    z = Interval(1-maxabs, maxabs-1);
  }

  //std::cout << "zero: [" << negpos.min << "," << negpos.max << "]\n";
  
  I = BiInterval(negpos.get_union(posneg), pospos.get_union(negneg), z);
 
  //std::cout << " ==> " << I << std::endl << std::endl;
  
  return I;
}


// X/Y = Z
// Y = X/Z
Mistral::BiInterval Mistral::BiInterval::anti_div_Y(const int arg) {

#ifdef _DEBUG_INTERVALS
  if(_DEBUG_INTERVALS) {
    std::cout << "compute anti div : " << arg << " / [min?,max?] = ";
  display(std::cout);
  std::cout << std::endl;
  }
#endif 

  BiInterval I(-INFTY, INFTY);

  if(arg) {

    Interval pos, neg, zpos, zneg;
    
    if(zero) {
      if(arg<0) {
	zpos = Interval(arg+1, INFTY);
	zneg = Interval(-INFTY, -arg-1);
      } else {
	zpos = Interval(-arg+1, INFTY);
	zneg = Interval(-INFTY, arg-1);
      }
    }

    if(arg > 0) {
      pos = positive.anti_div_Y_pos(arg);
      neg = negative.anti_div_Y_pos(arg);
    } else {
      neg = positive.anti_div_Y_neg(arg);
      pos = negative.anti_div_Y_neg(arg);
    }
   
    I = BiInterval(neg.get_union(zneg), pos.get_union(zpos), false);

  }

  return I;
}

Mistral::BiInterval Mistral::BiInterval::anti_div_Y(const Mistral::BiInterval arg) {

#ifdef _DEBUG_INTERVALS
  if(_DEBUG_INTERVALS) {
    std::cout << "compute anti div : " << arg << " / [min?,max?] = ";
  display(std::cout);
  std::cout << std::endl;
  }
#endif 

  BiInterval I(-INFTY, INFTY);
  Interval pos, neg, zpos, zneg;

  if(!(arg==0)) {

    Interval pospos = positive.anti_div_Y(arg.positive);
    
    //std::cout << "pospos: [" << pospos.min << "," << pospos.max << "]\n";
    
    Interval posneg = positive.anti_div_Y(arg.negative);
    
    //std::cout << "posneg: [" << posneg.min << "," << posneg.max << "]\n";
    
    Interval negpos = negative.anti_div_Y(arg.positive);
    
    //std::cout << "negpos: [" << negpos.min << "," << negpos.max << "]\n";
    
    Interval negneg = negative.anti_div_Y(arg.negative);
    
    //std::cout << "negneg: [" << negneg.min << "," << negneg.max << "]\n";
    
    Interval z;
    
    if(zero) {
      
      // int amin = 
      // int amax = arg.get_max_abs();
      int minabs = arg.get_min_abs(); //std::min(amin, amax);
      
      zpos = Interval(minabs+1, INFTY); 
      zneg = Interval(-INFTY, -minabs-1);
      
    }
    
    //std::cout << "zero: " << zneg << " u " << zpos << "\n";
    

    Interval K0 = negpos.get_union(posneg);

    //Interval K1 = negpos.get_union(posneg).get_union(zneg);
    Interval K1 = K0.get_union(zneg);

    Interval K2 = pospos.get_union(negneg).get_union(zpos);

    //std::cout << K0 << std::endl;
    
    //std::cout << K1 << std::endl;


    I = BiInterval(K1, K2, false);
   
  }

  //std::cout << " ==> " << I << std::endl << std::endl;

  return I;
}

Mistral::BiInterval Mistral::BiInterval::divided_by(const BiInterval arg, const Variable target) {
  BiInterval I;
  if(!(arg == 0)){
    if(*this == 0) {
      I = 0;
    } else { 
    
      Interval pospos = positive.divided_by(arg.positive, target);
      
      //  std::cout << "pospos: [" << pospos.min << "," << pospos.max << "]\n";
      
      Interval negneg = negative.divided_by(arg.negative, target);
      
      //  std::cout << "negneg: [" << negneg.min << "," << negneg.max << "]\n";
      
      
      Interval posneg = positive.divided_by(arg.negative, target);
      
      //  std::cout << "posneg: [" << posneg.min << "," << posneg.max << "]\n";
      
      Interval negpos = negative.divided_by(arg.positive, target);
      
      //  std::cout << "negpos: [" << negpos.min << "," << negpos.max << "]\n";
      
      
      // if this can be 0, then I can be 0
      I = BiInterval(posneg.get_union(negpos), pospos.get_union(negneg), zero);
    }
  }

  //std::cout << " ==> " << I << std::endl << std::endl;

  return I;
}

Mistral::BiInterval Mistral::BiInterval::operator/(const int arg) {
  Interval pos;
  Interval neg;

  // display(std::cout) ;
  // std::cout << " / " << arg << " = ";

  if(arg > 0) {
    pos = positive / arg;
    neg = negative / arg;
  } else if(arg < 0) {
    neg = positive / arg;
    pos = negative / arg;
  }

  BiInterval I(neg, pos, zero);

  // I.display(std::cout);
  // std::cout << std::endl << neg << std::endl << pos << std::endl;

  return I;
}



Mistral::BiInterval Mistral::BiInterval::operator*(const Mistral::BiInterval arg) {
  Interval pospos = positive * arg.positive;

  //  std::cout << "pospos: [" << pospos.min << "," << pospos.max << "]\n";

  Interval negneg = negative * arg.negative;

  //  std::cout << "negneg: [" << negneg.min << "," << negneg.max << "]\n";


  Interval posneg = positive * arg.negative;

  //  std::cout << "posneg: [" << posneg.min << "," << posneg.max << "]\n";

  Interval negpos = negative * arg.positive;

  //  std::cout << "negpos: [" << negpos.min << "," << negpos.max << "]\n";

  BiInterval I(posneg.get_union(negpos), pospos.get_union(negneg), zero || arg.zero);

  return I;
}

// the interval I such that this*I = arg     (I = this/arg)
Mistral::BiInterval Mistral::BiInterval::anti_mul(const Mistral::BiInterval arg) {
  BiInterval I;
  if(zero && arg.zero) {
    // if this and arg can be 0, then I is not constrained
    I = BiInterval(-INFTY, INFTY);
  } else if(arg == 0) {
    I = 0;
  } else {
    Interval pospos = positive.anti_mul(arg.positive);

    //    std::cout << "pospos: [" << pospos.min << "," << pospos.max << "]\n";

    Interval negneg = negative.anti_mul(arg.negative);

    //    std::cout << "negneg: [" << negneg.min << "," << negneg.max << "]\n";


    Interval posneg = positive.anti_mul(arg.negative);

    ///    std::cout << "posneg: [" << posneg.min << "," << posneg.max << "]\n";

    Interval negpos = negative.anti_mul(arg.positive);

    //    std::cout << "negpos: [" << negpos.min << "," << negpos.max << "]\n";


    // if arg cannot be 0, then I cannot be 0
    I = BiInterval(posneg.get_union(negpos), pospos.get_union(negneg), zero);
  }

  //  std::cout << " ==> " << I << std::endl << std::endl;

  return I;
}
// the interval I such that this/arg = I
Mistral::BiInterval Mistral::BiInterval::operator/(const Mistral::BiInterval arg) {
  BiInterval I;
  if(!(arg == 0)){
    if(*this == 0) {
      I = 0;
    } else { 
    
      Interval pospos = positive/(arg.positive);
      
      //  std::cout << "pospos: [" << pospos.min << "," << pospos.max << "]\n";
      
      Interval negneg = negative/(arg.negative);
      
      //  std::cout << "negneg: [" << negneg.min << "," << negneg.max << "]\n";
      
      
      Interval posneg = positive/(arg.negative);
      
      //  std::cout << "posneg: [" << posneg.min << "," << posneg.max << "]\n";
      
      Interval negpos = negative/(arg.positive);
      
      //  std::cout << "negpos: [" << negpos.min << "," << negpos.max << "]\n";
      
      
      // if this can be 0, then I can be 0
      I = BiInterval(posneg.get_union(negpos), pospos.get_union(negneg), zero);
    }
  }

  //std::cout << " ==> " << I << std::endl << std::endl;

  return I;
}

bool Mistral::BiInterval::operator==(const int x) const {
  if(x) {
    if(x<0) return x == negative.min && x == negative.max;
    else return x == positive.min && x == positive.max;
  } 
  return zero && positive.empty() && negative.empty();
}

void Mistral::BiInterval::operator=(const int x) {
  if(x) {
    if(x<0) {
      negative.min = x;
      negative.max = x;
      positive.min = +INFTY;
      positive.max = -INFTY;
      zero = false;
    } else {
      positive.min = x;
      positive.max = x;
      negative.min = +INFTY;
      negative.max = -INFTY;
      zero = false;
    }
  } else {
    positive.min = +INFTY;
    positive.max = -INFTY;
    negative.min = +INFTY;
    negative.max = -INFTY;
    zero = true;
  }
}

// Mistral::BiInterval Mistral::BiInterval::operator%(const int mod) {

// }

// Mistral::BiInterval Mistral::BiInterval::anti_modulo(const int mod);


Mistral::IntervalList::IntervalList() : Vector<Interval>() {}

Mistral::IntervalList::~IntervalList() {}

//void union_with(IntervalList& I);
void Mistral::IntervalList::intersect_with(const IntervalList& with, IntervalList& into) const {
  unsigned int current_self = 0;
  unsigned int current_with = 0;

#ifdef _DEBUG_INTERVALS
  if(_DEBUG_INTERVALS) {
    std::cout << "intersect\n";
  }
#endif

  Interval I, J;

  while(current_self < size && current_with < with.size) {

#ifdef _DEBUG_INTERVALS
  if(_DEBUG_INTERVALS) {
    std::cout << "self:  [" ;
    for(unsigned int i=current_self; i<size; ++i)
      std::cout << " " << stack_[i] ;
    std::cout << "]\nwith: [";
    for(unsigned int i=current_with; i<with.size; ++i)
      std::cout << " " << with[i] ;
    std::cout << "]\ninto: " << into << std::endl;
  }
#endif

    
    // find the next interval with maximum min
    if(with[current_with].min > stack_[current_self].min) {
      // we work with 'with[current_with]'
      // find the first interval in 'self' ending after the current interval in 'with' 
      while(current_self < size && stack_[current_self].max < with[current_with].min) ++current_self;
      if(current_self<size) {
	I = stack_[current_self];
	J = with[current_with];
      } else {
	//I.max = -INFTY;
	I.min = INFTY;
      }

    } else {
      // we work with 'with[current_with]'
      // find the first interval in 'with' ending after the current interval in 'self' 
      while(current_with < with.size && with[current_with].max < stack_[current_self].min) ++current_with;
      if(current_with<with.size) {
	I = with[current_with];    
	J = stack_[current_self];
      } else {
	//I.max = -INFTY;
	I.min = INFTY;
      }
      
    }



    // now we have max(I) >= min(J)
    // is there an overlap?
    if(// !I.empty() &&
       I.min <= J.max) {
      
#ifdef _DEBUG_INTERVALS
      if(_DEBUG_INTERVALS) {
	std::cout << " -> intersection: " << I << " ^ " << J << std::endl;
      }
#endif
      
      // put the intersection in 'into'
      into.add(Interval(std::max(I.min, J.min), 
			std::min(I.max, J.max)));
    }

    // increment the pointer(s)
    if(with[current_with].max >= stack_[current_self].max) ++current_self;
    else 
      //if(with[current_with].max <= stack_[current_self].max) 
      ++current_with;
  }
  
#ifdef _DEBUG_INTERVALS
  if(_DEBUG_INTERVALS) {
    std::cout << "self:  [" ;
    for(unsigned int i=current_self; i<size; ++i)
      std::cout << " " << stack_[i] ;
    std::cout << "]\nwith: [";
    for(unsigned int i=current_with; i<with.size; ++i)
      std::cout << " " << with[i] ;
    std::cout << "]\ninto: " << into << std::endl << std::endl;
  }
#endif
}


void Mistral::IntervalList::union_with(const IntervalList& with, IntervalList& into) const {
  unsigned int current_self = 0;
  unsigned int current_with = 0;

#ifdef _DEBUG_INTERVALS
  if(_DEBUG_INTERVALS) {
  std::cout << "union\n";
  }
#endif

  while(current_self < size || current_with < with.size) {

#ifdef _DEBUG_INTERVALS
  if(_DEBUG_INTERVALS) {
    std::cout << "self:  [" ;
    for(unsigned int i=current_self; i<size; ++i)
      std::cout << " " << stack_[i] ;
    std::cout << "]\nwith: [";
    for(unsigned int i=current_with; i<with.size; ++i)
      std::cout << " " << with[i] ;
    std::cout << "]\ninto: " << into << std::endl;
  }
#endif

    // which list has the next interval?
    int next_self = (current_self >= size ? INFTY : stack_[current_self].min);
    int next_with = (current_with >= with.size ? INFTY : with[current_with].min);


    // std::cout << "next_self = " << next_self << std::endl;
    // std::cout << "next_with = " << next_with << std::endl;

    
    Interval I;
    if(next_self < next_with) {
      I.min = next_self;
      I.max = stack_[current_self].max;
      ++current_self;
    } else {
      I.min = next_with;
      I.max = with[current_with].max;
      ++current_with;
    }

#ifdef _DEBUG_INTERVALS
  if(_DEBUG_INTERVALS) {
    std::cout << "  init I with " << I << std::endl;
  }
#endif

    // now find the end
    bool stop;
    do {
      stop = true;
      if(current_self < size && stack_[current_self].min <= I.max) {
	
#ifdef _DEBUG_INTERVALS
  if(_DEBUG_INTERVALS) {
	std::cout << "  overlap with " << stack_[current_self] << std::endl;
  }
#endif

	if(stack_[current_self].max > I.max) I.max = stack_[current_self].max;
	++current_self;
	stop = false;
      }
      if(current_with < with.size && with[current_with].min <= I.max) {

#ifdef _DEBUG_INTERVALS
  if(_DEBUG_INTERVALS) {
	std::cout << "  overlap with " << with[current_with] << std::endl;
  }
#endif

	if(with[current_with].max > I.max) I.max = with[current_with].max;
	++current_with;
	stop = false;
      }

#ifdef _DEBUG_INTERVALS
  if(_DEBUG_INTERVALS) {
      std::cout << "  -> " << I << std::endl;
  }
#endif

    } while(!stop);

    into.add(I);


    // if(next_self < next_with) {
    //   // the next start of an interval is in self, where is the end?

    //   // go through intervals in with, until the first that is strictly after current_self
    //   while(current_with < with.size && with[current_with].min <= stack_[current_self].max) ++current_with;
    //   into.add(Interval(next_self, (current_with ? std::max(stack_[current_self].max, with[current_with-1].max) : stack_[current_self].max)));
    //   ++current_self;
    
    // } else {
    //   // the next start of an interval is in self, where is the end?

    //   // go through intervals in with, until the first that is strictly after current_self
    //   while(current_self < size && stack_[current_self].min <= with[current_with].max) ++current_self;
    //   into.add(Interval(next_with, (current_self ? std::max(stack_[current_with].max, with[current_self-1].max) : with[current_with].max)));
    //   ++current_with;
    
    // }
  }

#ifdef _DEBUG_INTERVALS
  if(_DEBUG_INTERVALS) {
  std::cout << "self:  [" ;
  for(unsigned int i=current_self; i<size; ++i)
    std::cout << " " << stack_[i] ;
  std::cout << "]\nwith: [";
  for(unsigned int i=current_with; i<with.size; ++i)
    std::cout << " " << with[i] ;
  std::cout << "]\ninto: " << into << std::endl << std::endl;
  }
#endif
}


void Mistral::IntervalList::operator=(const IntervalList& l) {
  clear();
  for(unsigned int i=0; i<l.size; ++i) add(l[i]);
}


void Mistral::IntervalList::push(const int lb, const int ub) {
// #ifdef _CHECK_INTERVALS
//   if()
// #endif


  if(size && stack_[size-1].max >= lb-1)
    stack_[size-1].max = ub;
  else {
    Interval I(lb, ub);
    add(I);
  }
}


void Mistral::IntervalList::push(const Interval& I) {
// #ifdef _CHECK_INTERVALS
//   if()
// #endif


  if(size && stack_[size-1].max >= I.min-1)
    stack_[size-1].max = I.max;
  else {
    add(I);
  }
}

