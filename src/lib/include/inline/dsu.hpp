//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//////////////////////////////////////////////////////////////////////////////
namespace dsu {

inline char comp(char bp) {
   switch(bp) {
      case 'A' : return 'T'; 
      case 'C' : return 'G'; 
      case 'G' : return 'C'; 
      case 'T' : return 'A'; 
      case 'U' : return 'A'; 
      case 'Y' : return 'R'; 
      case 'R' : return 'Y'; 
      case 'K' : return 'M'; 
      case 'M' : return 'K'; 
      case 'B' : return 'V'; 
      case 'D' : return 'H'; 
      case 'H' : return 'D'; 
      case 'V' : return 'B'; 
      case 'N' : return 'N'; 
   case 'W' : return 'W';
   case 'S' : return 'S';
      case 'a' : return 't'; 
      case 'c' : return 'g'; 
      case 'g' : return 'c'; 
      case 't' : return 'a'; 
      case 'u' : return 'a'; 
      case 'y' : return 'r'; 
      case 'r' : return 'y'; 
      case 'k' : return 'm'; 
      case 'm' : return 'k'; 
      case 'b' : return 'v'; 
      case 'd' : return 'h'; 
      case 'h' : return 'd'; 
      case 'v' : return 'b'; 
      case 'n' : return 'n'; 
   case 'w' : return 'w';
   case 's' : return 's';
   	case '\0': return '\0';
   	default:
      	std::cerr<<"Unrecognized character: ["<<bp<<"] (dsu::comp)"<<std::endl;
      	return bp;
   }
}

}
