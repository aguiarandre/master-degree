#ifndef IGES_H

#define IGES_H

#include <stdio.h>
#include "native.h"
#include "BSpline.h"

#define FIELD_L 26  /*length of fixed length field in text file (sign, digits - double float) */


// Function Prototype:

void help(void);  
void format_number(char *n_s);                             
int read_iges_line(FILE *f, char *s);                /* read filerow from iges*/
void copy_field(char *s1, int poz, int x, char *s2); /* copy "x" chars from position "poz" from string "s1" to string "s2"*/
void parse_d_entry(char *ret1, char *ret2);          /* Parse single entry (pair of lines) in D section. */
void read_pd126_writetxt(char *s, int i, BSpline * spline, int splineId);
void read_pd128_writetxt(char *s, int begin);        /* read 128 PD in IGES and write to TXT */
int build_de_p(void);                                /* Build IGES DE section in memory - create points */
int build_de_g(void);                                /* Build IGES DE section in memory - create group of points */
int build_de_c(void);                                /* Build IGES DE section in memory - create ??? */
void writeiges_de(void);                             /* write DE section to IGES file */
void pd_linestring_complement(char *s, int ptr_de, int seq);
int is_numberlike(char znak);
void writeiges_pd(void);                             /* write PD section to IGES file */
void iges_read(int argc, char * igesFile, MeshStruct * mesh, BSpline * bSpline);

#endif //IGES_H
