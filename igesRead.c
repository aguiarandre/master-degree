#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "igesRead.h"
#include "BSpline.h"
#include "native.h"

#define debugging 0  /* 0 = normal mode, 1 = debugging */


FILE *fr, *fw, *fw2, *fr2, *fw3, *fw4;
/* options for i-mode */
int  o_write_3d,            /* 1 = write 3D coords */
     o_write_fixed_len,     /* 1 = write coords in fixed 26-char fields */
     o_write_dec_comma,     /* 1 = write decimal comma instead of dot */
     o_write_fixed_fp,      /* 1 = write fixed point format f10.4 (eliminate options "f" and ",") */
     o_write_float_fp,      /* 1 = write float point format g10.4 (eliminate options "f" and ",") */
     o_write_layername,     /* 1 = write name of the layer */
     o_extract_unblanked,   /* 1 = extract only unblanked */
     o_extract_pt,          /* 1 = extract points */
     o_extract_crv,         /* 1 = extract 126 curve */
     o_extract_surf;        /* 1 = extract 128 surface */
/* options for t-mode */
char o_crea_ent;            /* p, g, c for corresponding entity */
int  o_incr_layer;          /* 1 = increase layer*/

int entity[MAXENTITY][4];
  /* [][0] -entity_type(406-3=level,106,116=point,126=curve,128=surface)
   * [][1] -PDpointer
   * [][2] -PDcount (for export: PDcount= 106: entity_coords, 126: 1 +3*entity_coords +1; (1 header, n knot seq, n weights, n ctr pts, 1 unit normal)
   * [][3] -layer
   * */
   
enum entity_j_name{
  E_TYPE= 0,
  PD_PTR= 1,
  PD_CNT= 2,
  ELAYER= 3
};

// global vars:

int entity_sum=0,             /* A number of entities in array entity[] */
    emark=0;                  /* for cycling over entity[] */
int entity_coords[MAXENTITY]; /* A number of coordinates in entity */
#define MAX_LAYERNAME_LEN 79             /* maximum length of layer name string */
char layer[100][MAX_LAYERNAME_LEN];      /* max. 100 layers */
/* [x][ ] -layer Nr.
 * [ ][x] -layer name
 * */
int lmark=0;             			     /* for cycling over layer[] */



void help(void) {
  printf("IGSPREAD 0.5,**FREEWARE** Peter Gasparovic 2010-06-09\n");
  printf(" Description: \n");
  printf("  i-mode: program extracts coordinates of points (116, 106-2), and\n");
  printf("          control points of B-spline curves and surfaces (126, 128)\n");
  printf("          from IGES file and writes them to text file \"out.txt\"\n");
  printf("  t-mode: program reads coordinates in text file and creates IGES file\n");
  printf("          \"out.igs\". Coordinates must be blank separated, 2D or 3D,\n");
  printf("          each point on its own line. Empty line can increase layer.\n");
  printf("          Example of coords: -.454e02 +12,022E+3 \n");
  printf(" Synopsis:\n");
  printf("  igspread -i[3][f][,][u][p][c][s] <file.igs>\n");
  printf("             3 -write 3D coords (x,y,z), -default is 2D\n");
  printf("             f -write aligned to 26 char fields \n");
  printf("             , -write colon as decimal separator \n");
  printf("             x -write fixed point format 12.6f \n");
  printf("             g -write double point format 12.6g (eliminate option \"x\") \n");
  printf("             l -name of the layer before coordinates\n");
  printf("             u -extract only unblanked entities\n");
  printf("             p -add points (116,106) to output\n");
  printf("             c -add B-spline curves (126) to output, -default\n");
  printf("             s -add B-spline surfaces (128) to output\n");
  printf("  	       igspread -t[p|g|c][l] <file.txt>\n");
  printf("             p -create points (116) - default\n");
  printf("             g -create groups of points (106-2)\n");
  printf("             c -create curves (126-0, polyline)\n");
  printf("             l -increase layer after each entity\n\n");
}


void format_number(char *n_s) {
  /* format number string according options */
  /* char temporary[FIELD_L]; */
  int n_len;
  int i, j;
  double temp_fp;
  char temp_n_s[FIELD_L+1];

  n_len= strlen(n_s);
  /* regularize E,e,D,d -> e */
  for (i=0; i < n_len; i++) {
    if(n_s[i] == 'E' || n_s[i] == 'e' || n_s[i] == 'D' || n_s[i] == 'd') {
      n_s[i]='e';
    }
  }
  /* eliminate e0 e+0 e-0 e+000 at the end
   * cycle from the end - on 'e', move the end, else jump out
   * */
  if(n_len > 2) {
    if(n_s[n_len - 1] == '0') {
      for (i=n_len-1; i >= 0; i--) {
        switch (n_s[i]) {
          case '0':
          case '+':
          case '-':
            break;
          case 'e':
            n_s[i]= '\0';
            n_len= i;
            goto endfor_fn;
          default : goto endfor_fn;
        }
      }
      endfor_fn:  /* jump out of FOR from inside of SWITCH */
      ;
    }
  }

  /* convert string according different floating point format */
  if( o_write_fixed_fp == 1 || o_write_float_fp == 1 ) {
    temp_fp= atof(n_s);
    if( o_write_fixed_fp == 1) {
      sprintf(n_s, "%12.6f", temp_fp);
    }
    else if(o_write_float_fp == 1) {
      sprintf(n_s, "%12.6g", temp_fp);
    }
    n_len= strlen(n_s);
  }

  /* switch between decimal dot and comma */
  if ( o_write_dec_comma == 1 ) {
    for (i=0; i < n_len; i++) {
      if(n_s[i] == '.') {
        n_s[i]= ','; break;
      }
    }
  }

  /* complement from the left side the fixed field length */
  if (o_write_fixed_len == 1) {
    for (i=0; i<(FIELD_L -n_len); i++) {
      temp_n_s[i]=' ';
    }
    j=0;
    for (i=(FIELD_L -n_len); i<(FIELD_L); i++) {
      temp_n_s[i]=n_s[j];
      j++;
    }
    temp_n_s[FIELD_L]= '\0';
    strcpy(n_s, temp_n_s);
    n_len= strlen(n_s);
  }
}

int read_iges_line(FILE *f, char *s) {
  /* read filerow from iges*/
  int i;
  int c;
  for (i=0; i < 80; i++) {
    c=getc(f);
    if (c == EOF) return -1;
    s[i]=(char) c;
    if (i==0) { /* eat CR LF in line beginning and reset "i" to zero */
      switch (s[i]) {
        case '\x0D': i--; break;
        case '\x0A': i--; break;
      }
    }
  }
  return 0;
}

void copy_field(char *s1, int poz, int x, char *s2) {
  /* copy "x" chars from position "poz" from string "s1" to string "s2"*/
  int i, j;
  i=poz;
  j=0;
  while (j < x) {
    s2[j]=s1[i];
    i++;
    j++;
  }
}

void parse_d_entry(char *ret1, char *ret2) {
  /* Parse single entry (pair of lines) in D section. */
  char fld[9];                   /* for reading the value of array */
  int kod, form, pd_ptr, pd_count, layer;
  fld[8]='\0';
  /* Read:
   * entity type (A), form (B), pointer to PD line (C), PD count (D), layer (E)
   * and if the entity is recognised (406-3,106-2,116,126,128) load it into entity[][]
   * */
  copy_field(ret1,0,8,fld);  kod=     atoi(fld);       /* A */
  copy_field(ret2,32,8,fld); form=    atoi(fld);       /* B */
  copy_field(ret1,8,8,fld);  pd_ptr=  atoi(fld);       /* C */
  copy_field(ret2,24,8,fld); pd_count=atoi(fld);       /* D */
  copy_field(ret1,32,8,fld); layer=   atoi(fld);       /* E */
  switch (kod) {
    case 406:
      if (form != 3) break;                       /* only layers (406-3)*/
      entity[emark][E_TYPE]=kod;
      entity[emark][PD_PTR]=pd_ptr;
      entity[emark][PD_CNT]=pd_count;
      entity[emark][ELAYER]=layer;
      entity_sum++;
      emark++;
      break;
    case 106:
      if (o_extract_pt == 0) break;
      if (form != 2) break;                       /* only points (106-2)*/
      if (o_extract_unblanked == 1 && ret1[65] == '1') break; /* only unblanked ! */
      entity[emark][E_TYPE]=kod;
      entity[emark][PD_PTR]=pd_ptr;
      entity[emark][PD_CNT]=pd_count;
      entity[emark][ELAYER]=layer;
      entity_sum++;
      emark++;
      break;
    case 116:
      if (o_extract_pt == 0) break;
      if (o_extract_unblanked == 1 && ret1[65] == '1') break;
      entity[emark][E_TYPE]=kod;
      entity[emark][PD_PTR]=pd_ptr;
      entity[emark][PD_CNT]=pd_count;
      entity[emark][ELAYER]=layer;
      entity_sum++;
      emark++;
      break;
    case 126:
      if (o_extract_crv == 0) break;
      if (o_extract_unblanked == 1 && ret1[65] == '1') break;
      entity[emark][E_TYPE]=kod;
      entity[emark][PD_PTR]=pd_ptr;
      entity[emark][PD_CNT]=pd_count;
      entity[emark][ELAYER]=layer;
      entity_sum++;
      emark++;
      break;
    case 128:
      if (o_extract_surf == 0) break;
      if (o_extract_unblanked == 1 && ret1[65] == '1') break;
      entity[emark][E_TYPE]=kod;
      entity[emark][PD_PTR]=pd_ptr;
      entity[emark][PD_CNT]=pd_count;
      entity[emark][ELAYER]=layer;
      entity_sum++;
      emark++;
      break;
  }
}

// B-spline
void read_pd126_writetxt(char * s, int begin, BSpline * spline, int splineId)
{
  /* read 126 PD in IGES and write to TXT */

  static char x[FIELD_L+1], y[FIELD_L+1], z[FIELD_L+1], weight[FIELD_L+1], knot[FIELD_L+1];
  static enum {
    HEADER= 0,
    X= 1,
    Y= 2,
    Z= 3,
    TAIL_AFTER_LAST_PT= 4,
    W = 5,
    K = 6
  } phase;
  static int m;                  /* marker in strings x, y, z */
  static int pt_counter;         /* a number of ctrl. points */
  static int b_spline_param[2];    /* b[0]-maxindex_of_ctrl_pts, b[1]-degree_of_Bspl */
  enum bsp_i_name{
      CTRL_PT_MAX= 0,
      B_SPL_DEGREE= 1
    };
  char temp_string[5];
  int i, n;        /* i-order of the character, k-fix.write, n-temporary */
  static int j, jmax, kmax, wmax; /* j-order of head., jmax-length of header */
  if (begin == 1) {   /* first PD line of entity */
    phase=0;
    pt_counter=0; 
    m=0;
    j=0; jmax=3; kmax = 3; wmax =3 ;      /* jmax is definitely computed in the step j=2 */
    
    n=0;
  }
  for (i=0; i <= 64; i++) {
    switch (phase) {
      case HEADER:
        switch (j) {
          case 1:
            if (s[i] == ',') {
              temp_string[n]='\0';
              b_spline_param[CTRL_PT_MAX]=atoi(temp_string);
	      j++; n=0;
            }
            else {
              temp_string[n]=s[i]; n++;
            }
            break;
          case 2:
            if (s[i] == ',') {
              temp_string[n]='\0';
              b_spline_param[B_SPL_DEGREE]=atoi(temp_string);
              j++; n=0; // so j is the counter of commas (,)
	      
	      int kk, nn, aa, mm;
	      kk = b_spline_param[0];
	      mm = b_spline_param[1];
	      spline[ splineId ].k = kk;
          spline[ splineId ].m = mm;
          spline[ splineId ].bSplineId = splineId;
          //params_global[splineId][0] = b_spline_param[0];
	      //params_global[splineId][1] = b_spline_param[1];
	      nn = 1 + kk - mm;
	      aa = nn + 2*mm;

	      kmax = 7+aa;
	      wmax = 8+aa+kk;
              jmax=b_spline_param[CTRL_PT_MAX]+b_spline_param[CTRL_PT_MAX]+b_spline_param[B_SPL_DEGREE]+10; // My idea is to develop some kind of mechanism closely related to this jmax.
            }
            else {
              temp_string[n]=s[i]; n++;
            }
            break;
          default:
            if (s[i] == ',') j++; 
	    if (j >= 7 && j <= kmax) {
phase = K;
fprintf(fw,"-------------------- B-Spline # %d --------------------\n \nKnot sequence:\n", splineId+1);  /* This is just to make it easier to visualize - not practical at all if I want to read the file afterwards, but ... */
} 

	}
        break;
/* Same thing done to compute the knot vector, weights and control points. They are different 'phase' cases, though */

// Knot Vector 

      case K:
        if (s[i] == ' ')                  
          goto endline126;
        else if (s[i] == ',') 
        {                 
          knot[m]= '\0';
          format_number(knot);/** TODO */
      spline[ splineId ].knots[ j - 7 ] = atof(knot); 
      //knots[bspline_counter][j-7] = atof(knot); // when you reach a comma, you have a full number written in string format. Now I can translate it to a float. (doubts - precision?)
	  
      fprintf(fw, "%s\n", knot);
          m= 0;
	  j++; 			/* j counts the number of numbers read so far (j = 0 -> '126') */
	  if (j == kmax+1) { //if it reach the position of the weights,
phase = W;
fprintf(fw,"\nWeights\n");
} 
          break;
        }
        else {                                  
          knot[m]= s[i];
	  
          m++;
          break;
        }
        break;

// Weights

      case W:
      if (s[i] == ' ')                     
        goto endline126;
    
      else if (s[i] == ',') 
      {                 
          weight[m]= '\0';
          format_number(weight);
	  spline[ splineId ].weights[ j - 8 - (1+b_spline_param[0]+b_spline_param[1])] = atof( weight ); 
      //weights[bspline_counter][j-8-(1+b_spline_param[0]+b_spline_param[1])] = atof(weight); //global variable.
	  fprintf(fw, "%s\n", weight);
          m= 0;
	  j++;
	  if (j == jmax) { //if it reach the position of the control points,
phase = X; 
fprintf(fw,"\nControl Points\n");
}
          break;
        }
        else {                                  
          weight[m]= s[i];
	  
          m++;
          break;
        }
        break;
	    
// Control Points.
      case X:
        if (s[i] == ' ') {                      /* end of PD line */
          goto endline126;
        }
        else if (s[i] == ',') {                 /* next coordinate */
          x[m]= '\0';
          format_number(x);
          spline[ splineId ].x_ctrl[j-jmax] = atof(x);
	  //x_global[bspline_counter][j-jmax] = atof(x);
          m= 0;
          phase= Y;
          break;
        }
        else {                                  
          x[m]= s[i];
          m++;
          break;
        }
        break;
      case Y:
        if (s[i] == ' ') {
          goto endline126;
        }
        else if (s[i] == ',') {
          y[m]= '\0';
          format_number(y);
      spline[ splineId ].y_ctrl[j-jmax] = atof(y);
	  //y_global[bspline_counter][j-jmax] = atof(y);
          m= 0;
          phase= Z;
          break;
        }
        else {
          y[m]= s[i];
          m++;
          break;
        }
        break;
      case Z:
        if (s[i] == ' ') {                      /* end of PD line */
          goto endline126;
        }
        else if (s[i] == ',') {
          z[m]='\0';
          format_number(z);
      spline[ splineId ].z_ctrl[j-jmax] = atof(z);
	  //z_global[bspline_counter][j-jmax] = atof(z);
          fprintf(fw, "%s%s%s\n", x, y, z);
          pt_counter++;
          if (pt_counter < (b_spline_param[CTRL_PT_MAX]+1)) {
            m= 0;
            phase=X;
	  j++;
          }
          else {                                               /* Z of the last point */
            entity_coords[emark]=pt_counter;
            phase=TAIL_AFTER_LAST_PT;
          }
          break;
        }
        else {                                
          z[m]= s[i];
          m++;
          break;
        }
        break;
      case TAIL_AFTER_LAST_PT:
        goto endline126;
        break;
    }
  }
  endline126:
  ;
}

void read_pd128_writetxt(char *s, int begin) {
  /* read 128 PD in IGES and write to TXT */
  static char x[FIELD_L+1], y[FIELD_L+1], z[FIELD_L+1];
  static int phase; /* 0-header, 1-x, 2-y, 3-z */
  static int m;       /* marker in strings x, y, z */
  static int c;    /* counter of points */
  static int b[4]; /* b[0,1]-maxindex_of_ctrl_pts, b[2,3]-degree_of_Bspl */
  char b_s[5];
  int i;  /* i-order of the character, k-fix.write,  */
  static int n, j, jmax;  /* n-temporary, j-order of head., jmax-length of header */
  if (begin == 1) {
    phase=0;
    c=0;
    m=0;
    j=0; jmax=5;  /* jmax is definitely computed only when j=4 */
    n=0;
  }
  for (i=0; i <= 64; i++) {
    switch (phase) {
      case 0:
        switch (j) {
          case 1:
          case 2:
          case 3:
            if (s[i] == ',') {
              b_s[n]='\0';
              b[j-1]=atoi(b_s);
              j++; n=0;
            }
            else {
              b_s[n]=s[i]; n++;
            }
            break;
          case 4:
            if (s[i] == ',') {
              b_s[n]='\0';
              b[j-1]=atoi(b_s);
              j++; n=0;
              jmax=b[0]+b[1]+b[2]+b[3]+14+(b[0]+1)*(b[1]+1);
            }
            else {
              b_s[n]=s[i]; n++;
            }
            break;
          default:
            if (s[i] == ',') j++;
            if (j == jmax) phase=1;
        }
        break;
      case 1:
        if (s[i] == ' ') {                      /* end of PD line */
          goto endline128;
        }
        else if (s[i] == ',') {                 /* next coordinate */
          x[m]= '\0';
          format_number(x);
          m= 0;
          phase= 2;
          break;
        }
        else {                                  /* digit + - E e D d */
          x[m]= s[i];
          m++;
          break;
        }
        break;
      case 2:
        if (s[i] == ' ') {
          goto endline128;
        }
        else if (s[i] == ',') {
          y[m]= '\0';
          format_number(y);
          m= 0;
          phase= 3;
          break;
        }
        else {
          y[m]= s[i];
          m++;
          break;
        }
        break;
      case 3:
        if (s[i] == ' ') {                      /* end of PD line */
          goto endline128;
        }
        else if (s[i] == ',' || s[i] == ';') {  /* the last coordinate */
          z[m]='\0';
          format_number(z);
          if (o_write_3d == 0) fprintf(fw, "%s %s\n", x, y);
          else fprintf(fw, "%s %s %s\n", x, y, z);
          m= 0;
          c++;
          if (c == (b[0]+1)*(b[1]+1)) phase=4; else phase=1;
          if (c % (b[0]+1) == 0) fprintf(fw, "\n");
          break;
        }
        else {                                  /* digit + - E e D d */
          z[m]= s[i];
          m++;
          break;
        }
        break;
      case 4:
        goto endline128;
    }
  }
  endline128:
  ;
}

int build_de_p(void) {
  /* Build IGES DE section in memory - create points - 116 */
  int  current_line_chars,    /* to detect non-empty line (number of non-blank characters in line) */
       previous_line_chars,   /* to detect that the current line is not the first (number of non-blank characters in previous line) */
       param_seq;             /* pointer to PD data */
  char character;
  static int lay;
  current_line_chars=0;
  previous_line_chars=0;
  param_seq=1;
  lay=0;
  character=' ';
  fseek(fr, 0L, SEEK_SET);
  while(character != EOF) {
    character=getc(fr);
    switch (character) {
      case '\n':
        if (current_line_chars != 0) {
          entity_coords[emark]=1;
          entity[emark][E_TYPE]=116;
          entity[emark][PD_PTR]=param_seq;
          entity[emark][PD_CNT]=1;
          entity[emark][ELAYER]=lay;
          param_seq+=entity[emark][PD_CNT];
          entity_sum++;
          if (entity_sum == MAXENTITY) return -1;
          emark++;
        }
        else {
          if (previous_line_chars != 0) {
            if (o_incr_layer == 1) lay++;
          }
        }
        previous_line_chars=current_line_chars;
        current_line_chars=0;
        break;
      case ' ':
      case '\t':
      case '\r': break;
      default: current_line_chars++; break;
    }
  }
  return 0;
}

int build_de_g(void) {
  /* Build IGES DE section in memory - create group of points - 106 */
  int  current_line_chars,    /* to detect non-empty line (number of non-blank characters in line) */
       previous_line_chars,   /* to detect that the current line is not the first (number of non-blank characters in previous line) */
       coordinates,           /* number of coordinates (contiguous lines) */
       param_seq;             /* pointer to PD data */
  char character;
  static int lay;
  current_line_chars=0;
  previous_line_chars=0;
  coordinates=0;
  param_seq=1;
  lay=0;
  character=' ';
  fseek(fr, 0L, SEEK_SET);
  do {
    character=getc(fr);
    switch (character) {
      case '\n':
        if (current_line_chars != 0) {                       /* \n after non-empty line (with coordinates) */
          coordinates++;
        }
        else if (previous_line_chars != 0) {                 /* \n in empty line after contiguous block */
          entity_coords[emark]=coordinates;
          entity[emark][E_TYPE]=106;
          entity[emark][PD_PTR]=param_seq;
          entity[emark][PD_CNT]=coordinates;
          entity[emark][ELAYER]=lay;
          param_seq+=entity[emark][PD_CNT];
          coordinates=0;
          entity_sum++;
          if (entity_sum == MAXENTITY) return -1;
          emark++;
          if (o_incr_layer == 1) lay++;
        }
        previous_line_chars=current_line_chars;
        current_line_chars=0;
        break;
      case EOF:
        if (previous_line_chars != 0) {
          entity_coords[emark]=coordinates;
          entity[emark][E_TYPE]=106;
          entity[emark][PD_PTR]=param_seq;
          entity[emark][PD_CNT]=coordinates;
          entity[emark][ELAYER]=lay;
          param_seq+=entity[emark][PD_CNT];
          coordinates=0;
          entity_sum++;
          if (entity_sum == MAXENTITY) return -1;
          emark++;
          if (o_incr_layer == 1) lay++;
        }
        break;
      case ' ':
      case '\t':
      case '\r': break;
      default: current_line_chars++; break;            /* detect non-empty line */
    }
  } while (character != EOF);
  return 0;
}

int build_de_c(void) {
  /* Build IGES DE section in memory - create curve - 126*/
  int  current_line_chars,    /* to detect non-empty line (number of non-blank characters in line) */
  previous_line_chars,   /* to detect that the current line is not the first (number of non-blank characters in previous line) */
  coordinates,           /* number of coordinates (contiguous lines) */
  param_seq;             /* pointer to PD data */
  char character;
  static int lay;
  current_line_chars=0;
  previous_line_chars=0;
  coordinates=0;
  param_seq=1;
  lay=0;
  character=' ';
  fseek(fr, 0L, SEEK_SET);
  do {
    character=getc(fr);
    switch (character) {
      case '\n':
        if (current_line_chars != 0) {                       /* \n after non-empty line (with coordinates) */
          coordinates++;
        }
        else if (previous_line_chars != 0) {                 /* \n in empty line after contiguous block */
          entity_coords[emark]=coordinates;
          entity[emark][E_TYPE]=126;
          entity[emark][PD_PTR]=param_seq;
          entity[emark][PD_CNT]=2+3*coordinates;
          entity[emark][ELAYER]=lay;
          param_seq+=entity[emark][PD_CNT];
          coordinates=0;
          entity_sum++;
          if (entity_sum == MAXENTITY) return -1;
          emark++;
          if (o_incr_layer == 1) lay++;
        }
        previous_line_chars=current_line_chars;
        current_line_chars=0;
        break;
      case EOF:
        if (previous_line_chars != 0) {
          entity_coords[emark]=coordinates;
          entity[emark][E_TYPE]=126;
          entity[emark][PD_PTR]=param_seq;
          entity[emark][PD_CNT]=2+3*coordinates;
          entity[emark][ELAYER]=lay;
          param_seq+=entity[emark][PD_CNT];
          coordinates=0;
          entity_sum++;
          if (entity_sum == MAXENTITY) return -1;
          emark++;
          if (o_incr_layer == 1) lay++;
        }
        break;
      case ' ':
      case '\t':
      case '\r': break;
      default: current_line_chars++; break;            /* detect non-empty line */
    }
  } while (character != EOF);
  return 0;
}

void pd_linestring_complement(char *s, int ptr_de, int seq) {
    
  int i;
  char endstring[20];
  for (i=strlen(s); i <= 64; i++) {
    s[i]=' ';
  }
  s[65]='\0';
  sprintf(endstring, "%07dP%7d", ptr_de, seq);
  strcat(s, endstring);
}

int is_numberlike(char znak) {
  if ( (znak >= '0' && znak <= '9')
      || znak == '-'
        || znak == '+'
          || znak == '.'
            || znak == ','
              || znak == 'E'
                || znak == 'e'
                  || znak == 'D'
                    || znak == 'd' ) {
    return 1;
  }
  else {
    return 0;
  }
}

void writeiges_pd(void) {
    
  /* write PD section to IGES file */
  int i;
  char znak; // tradução: Mark
  int  i_znak;

  char riadok[81]; //tradução: Line
  char zbytok[81]; //tradução: Residue

  char x[FIELD_L+1], y[FIELD_L+1], z[FIELD_L+1];
  char *xyz;      /* pointer to x[], y[] or z[] */
  int seq;        /* PD line position (seq number) */
  int line_count; /* number of coordinates for current entity */

  /* State variable for 116 and 106 */
  int complet;

  /* State variables for 126 */
  int prev_char_numberlike;               /* 0,1  Is the previous character number-like? */
  int curr_char_numberlike;               /* 0,1  Is the current character number-like? */
  int numberlike_in_line_counter;         /* 0,1,2,... How many characters in line are number-like? */
  enum {NOTHING, X, Y, Z} the_last_coord; /* last coord being started (0= no coordinate yet), 1 = X, 2 = Y, 3 = Z */
  int prev_line_nonblank;                 /* 0,1  Contains the previous line number-like characters? */
  int curr_line_nonblank;                 /* 0,1  Contains the current line number-like characters? */

  /* Initializers - global */
  emark= 0;
  fseek(fr, 0L, SEEK_SET);

  /* Initializers - local */
  znak= ' ';
  i_znak= 0;
  riadok[0]= '\0';
  xyz= x;
  seq= 1;
  line_count=1;

  complet=0;

  prev_char_numberlike= 0;
  curr_char_numberlike= 0;
  numberlike_in_line_counter= 0;
  the_last_coord= NOTHING;
  prev_line_nonblank= 0;
  curr_line_nonblank= 0;

//------------------------------------------------------------------------------POINTS 116-----------------------------------------------------------------//
  switch (o_crea_ent) {
    case 'p': {                 /* pts - 116 */
      while (znak != EOF) {
        znak=getc(fr);
        switch (complet) {
          case 0:
            if ((znak >= '0' && znak <= '9') || znak == '-' || znak == '+' ||
                znak == '.') {
              sprintf(riadok, "%d,", 116);
              complet=1;
              x[i_znak]=znak;
              i_znak++;
            }
            else {
              if (znak == ',') {
                sprintf(riadok, "%d,", 116);
                complet=1;
                x[i_znak]='.';
                i_znak++;
              }
            }
            break;
          case 1:
            if ((znak >= '0' && znak <= '9') || znak == '-' || znak == '+' ||
                znak == '.') {
              x[i_znak]=znak;
              i_znak++;
            }
            else {
              switch (znak) {
                case 'e':
                case 'E':
                  x[i_znak]='D';
                  i_znak++;
                  break;
                case ',':
                  x[i_znak]='.';
                  i_znak++;
                  break;
                case ' ':
                case '\t':
                  x[i_znak]='\0';
                  if (i_znak != 0) {
                    complet=2;
                    i_znak=0;
                  }
                  break;
                case '\n':
                  x[i_znak]='\0';
                  strcpy(y, "0D0");
                  strcpy(z, "0D0");
                  sprintf(zbytok, "%s,%s,%s,0;", x, y, z);
                  strcat(riadok, zbytok);
                  pd_linestring_complement(riadok, 2*emark+1, seq);
                  fprintf(fw, "%s\n", riadok);
                  seq++;
                  line_count++;/* condition for complex entities */
                  emark++;
                  complet=0;
                  i_znak=0;
                  break;
              }
            }
            break;
          case 2:
            if ((znak >= '0' && znak <= '9') || znak == '-' || znak == '+' ||
                znak == '.') {
              y[i_znak]=znak;
              i_znak++;
            }
            else {
              switch (znak) {
                case 'e':
                case 'E':
                  y[i_znak]='D';
                  i_znak++;
                  break;
                case ',':
                  y[i_znak]='.';
                  i_znak++;
                  break;
                case ' ':
                case '\t':
                  y[i_znak]='\0';
                  if (i_znak != 0) {
                    complet=3;
                    i_znak=0;
                  }
                  break;
                case '\n':
                  y[i_znak]='\0';
                  if (i_znak == 0) strcpy(y, "0D0");
                  strcpy(z, "0D0");
                  sprintf(zbytok, "%s,%s,%s,0;", x, y, z);
                  strcat(riadok, zbytok);
                  pd_linestring_complement(riadok, 2*emark+1, seq);
                  fprintf(fw, "%s\n", riadok);
                  seq++;
                  line_count++;/* condition for complex entities */
                  emark++;
                  complet=0;
                  i_znak=0;
                  break;
              }
            }
            break;
          case 3:
            if ((znak >= '0' && znak <= '9') || znak == '-' || znak == '+' ||
                znak == '.') {
              z[i_znak]=znak;
              i_znak++;
            }
            else {
              switch (znak) {
                case 'e':
                case 'E':
                  z[i_znak]='D';
                  i_znak++;
                  break;
                case ',':
                  z[i_znak]='.';
                  i_znak++;
                  break;
                case ' ':
                case '\t':
                  z[i_znak]='\0';
                  break;
                case '\n':
                  z[i_znak]='\0';
                  if (i_znak == 0) strcpy(z, "0D0");
                  sprintf(zbytok, "%s,%s,%s,0;", x, y, z);
                  strcat(riadok, zbytok);
                  pd_linestring_complement(riadok, 2*emark+1, seq);
                  fprintf(fw, "%s\n", riadok);
                  seq++;
                  line_count++;/* condition for complex entities */
                  emark++;
                  complet=0;
                  i_znak=0;
                  break;
              }
            }
            break;
        }
      }
      break;
    } /*end case p*/

//-------------------------------------------------------------------------GROUPS OF POINTS 106--------------------------------------------------------------//
    case 'g': {                 /* groups of pts - 106 */
      while (znak != EOF) {
        znak=getc(fr);
        switch (complet) {
          case 0:
            if ((znak >= '0' && znak <= '9') || znak == '-' || znak == '+' ||
                 znak == '.') {
              sprintf(riadok, "106,2,%d,", entity[emark][PD_CNT]);
              complet=1;
              x[i_znak]=znak;
              i_znak++;
            }
            else {
              if (znak == ',') {
                sprintf(riadok, "106,2,%d,", entity[emark][PD_CNT]);
                complet=1;
                x[i_znak]='.';
                i_znak++;
              }
            }
            break;
          case 1:
            if ((znak >= '0' && znak <= '9') || znak == '-' || znak == '+' ||
                 znak == '.') {
              x[i_znak]=znak;
              i_znak++;
            }
            else {
              switch (znak) {
                case 'e':
                case 'E':
                  x[i_znak]='D';
                  i_znak++;
                  break;
                case ',':
                  x[i_znak]='.';
                  i_znak++;
                  break;
                case ' ':
                case '\t':
                  x[i_znak]='\0';
                  if (i_znak != 0) {
                    complet=2;
                    i_znak=0;
                  }
                  break;
                case '\n':
                  x[i_znak]='\0';
                  strcpy(y, "0D0");
                  strcpy(z, "0D0");
                  if (line_count == entity[emark][PD_CNT]) {
                    sprintf(zbytok, "%s,%s,%s;", x, y, z);
                    strcat(riadok, zbytok);
                    pd_linestring_complement(riadok, 2*emark+1, seq);
                    fprintf(fw, "%s\n", riadok);
                    complet=0;
                    line_count=1;
                    emark++;
                  }
                  else {
                    sprintf(zbytok, "%s,%s,%s,", x, y, z);
                    strcat(riadok, zbytok);
                    pd_linestring_complement(riadok, 2*emark+1, seq);
                    fprintf(fw, "%s\n", riadok);
                    riadok[0]='\0';
                    line_count++;
                  }
                  seq++;
                  i_znak=0;
                  break;
              }
            }
            break;
          case 2:
            if ((znak >= '0' && znak <= '9') || znak == '-' || znak == '+' ||
                 znak == '.') {
              y[i_znak]=znak;
              i_znak++;
            }
            else {
              switch (znak) {
                case 'e':
                case 'E':
                  y[i_znak]='D';
                  i_znak++;
                  break;
                case ',':
                  y[i_znak]='.';
                  i_znak++;
                  break;
                case ' ':
                case '\t':
                  y[i_znak]='\0';
                  if (i_znak != 0) {
                    complet=3;
                    i_znak=0;
                  }
                  break;
                case '\n':
                  y[i_znak]='\0';
                  if (i_znak == 0) strcpy(y, "0D0");
                  strcpy(z, "0D0");
                  if (line_count == entity[emark][PD_CNT]) {
                    sprintf(zbytok, "%s,%s,%s;", x, y, z);
                    strcat(riadok, zbytok);
                    pd_linestring_complement(riadok, 2*emark+1, seq);
                    fprintf(fw, "%s\n", riadok);
                    complet=0;
                    line_count=1;
                    emark++;
                  }
                  else {
                    sprintf(zbytok, "%s,%s,%s,", x, y, z);
                    strcat(riadok, zbytok);
                    pd_linestring_complement(riadok, 2*emark+1, seq);
                    fprintf(fw, "%s\n", riadok);
                    complet=1;
                    riadok[0]='\0';
                    line_count++;
                  }
                  seq++;
                  i_znak=0;
                  break;
              }
            }
            break;
          case 3:
            if ((znak >= '0' && znak <= '9') || znak == '-' || znak == '+' ||
                 znak == '.') {
              z[i_znak]=znak;
              i_znak++;
            }
            else {
              switch (znak) {
                case 'e':
                case 'E':
                  z[i_znak]='D';
                  i_znak++;
                  break;
                case ',':
                  z[i_znak]='.';
                  i_znak++;
                  break;
                case ' ':
                case '\t':
                  z[i_znak]='\0';
                  break;
                case '\n':
                  z[i_znak]='\0';
                  if (i_znak == 0) strcpy(z, "0D0");
                  if (line_count == entity[emark][PD_CNT]) {
                    sprintf(zbytok, "%s,%s,%s;", x, y, z);
                    strcat(riadok, zbytok);
                    pd_linestring_complement(riadok, 2*emark+1, seq);
                    fprintf(fw, "%s\n", riadok);
                    complet=0;
                    line_count=1;
                    emark++;
                  }
                  else {
                    sprintf(zbytok, "%s,%s,%s,", x, y, z);
                    strcat(riadok, zbytok);
                    pd_linestring_complement(riadok, 2*emark+1, seq);
                    fprintf(fw, "%s\n", riadok);
                    riadok[0]='\0';
                    complet=1;
                    line_count++;
                  }
                  seq++;
                  i_znak=0;
                  break;
              }
            }
            break;
        }
      }
      break;
    } /*end case g*/

//---------------------------------------------------------------- CURVES - 126 -------------------------------------------------------------//

    case 'c': {                 /* curves - 126 */
      /* State variables: -global */
      emark= 0;
      fseek(fr, 0L, SEEK_SET);

      /* State variables: -local */
      znak= ' ';
      i_znak= 0;
      riadok[0]= '\0';
      xyz= x;
      seq= 1;
      line_count= 0;                  /* ! - different than for 116 and 106*/
      prev_char_numberlike= 0;
      curr_char_numberlike= 0;
      numberlike_in_line_counter= 0;
      the_last_coord= NOTHING;        /* {NOTHING, X, Y, Z} */
      prev_line_nonblank= 0;
      curr_line_nonblank= 0;


      while (znak != EOF) { // until the end of the file
        znak=getc(fr);      // znak reads the characters of the file


        /* At first, proces chars */
        prev_char_numberlike= curr_char_numberlike; 
        curr_char_numberlike= is_numberlike(znak); //check to see if this character is a number. (+ - , . e E d D 0 9)


        if (curr_char_numberlike == 1) { // if it is a number,
          numberlike_in_line_counter++;  /* count how many valid chars */


          /* replace non-standard characters */
          switch (znak) {
            case 'e':
            case 'E':
            case 'd':
            case 'D':
              znak= 'D'; //exchange e E d for D
              break;
            case ',':
              znak= '.'; //exchange , for .
              break;
          }
        }

        if (prev_char_numberlike == 0 && curr_char_numberlike == 1) {         /* chars pattern 0 1 */
          the_last_coord++;
          i_znak= 0;
          switch (the_last_coord) {
            case X:
              xyz=x;
              break;
            case Y:
              xyz=y;
              break;
            case Z:
              xyz=z;
              break;
            case NOTHING:
              break;
          }
          xyz[i_znak]=znak;
          i_znak++;
        }
        else if (prev_char_numberlike == 1 && curr_char_numberlike == 1) {    /* chars pattern 1 1 */
          xyz[i_znak]=znak;
          i_znak++;
        }
        else if (prev_char_numberlike == 1 && curr_char_numberlike == 0) {    /* chars pattern 1 0 */
          xyz[i_znak]='\0';
          i_znak++;
        }


        /* on \n and EOF */
        if (znak == '\n' || znak == EOF) {
          prev_line_nonblank = curr_line_nonblank;
          if (numberlike_in_line_counter > 0) {           /* Is current line non-blank? */
            curr_line_nonblank= 1;
            line_count++;
          }
          else {
            curr_line_nonblank= 0;
          }
          numberlike_in_line_counter= 0;                    /* already used */
          prev_char_numberlike = 0;                         /* ready for next line */
          curr_char_numberlike = 0;                         /* ready for next line */

          if (prev_line_nonblank == 0 && curr_line_nonblank == 1) {             /* write header */
            sprintf(riadok, "126,%d,1,0,0,1,0,", entity_coords[emark] -1);                   /* line 1 */
            pd_linestring_complement(riadok, 2*emark+1, seq); // riadok, número ímpar (um sim um não), 
            fprintf(fw, "%s\n", riadok);
            seq++;
            sprintf(riadok, "0.0D0,0.0D0,");                                              /* knot seq - the first line */
            pd_linestring_complement(riadok, 2*emark+1, seq);
            fprintf(fw, "%s\n", riadok);
            seq++;
            for (i=1; i < entity_coords[emark] -1; i++) {
              sprintf(riadok, "%d.0D0,", i);                                              /*      - intermediate lines */
              pd_linestring_complement(riadok, 2*emark+1, seq);
              fprintf(fw, "%s\n", riadok);
              seq++;
            }
            sprintf(riadok, "%d.0D0,%d.0D0,", entity_coords[emark] -1, entity_coords[emark] -1);  /*      - the last line */
            pd_linestring_complement(riadok, 2*emark+1, seq);
            fprintf(fw, "%s\n", riadok);
            seq++;
            for (i=1; i <= entity_coords[emark]; i++) {
              sprintf(riadok, "%d.D0,", i);                                         /* weights */
              pd_linestring_complement(riadok, 2*emark+1, seq);
              fprintf(fw, "%s\n", riadok);
              seq++;
            }
          }

          if (curr_line_nonblank == 1) {                                        /* write xyz */
            switch (the_last_coord) {
              case X:
                strcpy(y, "0D0");
                strcpy(z, "0D0");
                break;
              case Y:
                strcpy(z, "0D0");
                break;
              case Z:
                break;
              case NOTHING:
                break;
            }
            sprintf(riadok, "%s,%s,%s,", x, y, z);
            pd_linestring_complement(riadok, 2*emark+1, seq);
            fprintf(fw, "%s\n", riadok);
            seq++;
            the_last_coord= NOTHING;
            i_znak= 0;
          }

          if ((prev_line_nonblank == 1 && curr_line_nonblank == 0)
              || (znak == EOF && curr_line_nonblank == 1) ) {                   /* write tail */
            sprintf(riadok, "0.0D0,%d.0D0,0.0D0,0.0D0,0.0D0;", entity_coords[emark] -1 -1);                   /* line 1 */
            pd_linestring_complement(riadok, 2*emark+1, seq);
            fprintf(fw, "%s\n", riadok);
            seq++;
            emark++;
          }
        } /* end if (znak == '\n' || znak == EOF) */
      } /* end while (znak != EOF) */
      break;
    } /* end case c */
  } /* end switch (o_crea_ent) */
}

void iges_read(int argc, char * igesFile, MeshStruct* mesh, BSpline * bSpline)
{ 
    	
  int  i, count,k,
       err_code;
  
  char line[81],          /* text content of current line from IGES */
       pairline[81];      /* pair line (in D section of IGES file) */
  char pd_seq_s[8];
  int bspline_counter=0; /* count the number of 126-entities */
  int  pd_seq;           /* parameter data sequence number */
 // int (*build_de)(void); /* pointer to function writing the DE */
  int npts;

  double parameter1, parameter2;
  int index = 0, j;
  int splineId = 0;

  npts = 100; npts = npts - 1;


//------------------------------------------------------------------------------------------------------------------//
  
  line[80]=   '\0';
  pairline[80]='\0';
  pd_seq_s[7]=  '\0';
  
/* defaults in i-mode */
  o_write_3d=1;          /* 1 = write 3D coords */
  o_write_fixed_len=0;   /* 1 = write coords in fixed 26-char fields */
  o_write_dec_comma=0;   /* 1 = write decimal comma instead of dot */
  o_write_fixed_fp=1,    /* 1 = write fixed point format f10.4 (eliminate options "f" and ",") */
  o_write_float_fp=0,    /* 1 = write float point format g10.4 (eliminate options "f" and ",") */
  o_write_layername=0;   /* 1 = write name of the layer */
  o_extract_unblanked=0; /* 1 = extract only unblanked */
  o_extract_pt=0;        /* 1 = extract points */
  o_extract_crv=1;       /* 1 = extract 126 curve */
  o_extract_surf=0;      /* 1 = extract 128 surface */
/* defaults in t-mode */
  o_crea_ent='p';        /* p, g, c for corresponding entity */
  o_incr_layer=0;        /* 1 = increase layers */

/* Analyze command line arguments */
  if (argc == 1) {help(); goto end;}
  
  
/*--------------------- IGES mode -----------------------------------------*/
  // Read Iges File.
  fr= fopen( igesFile , "rb");
  if (fr == NULL) goto error3;
  fw= fopen("out.txt", "w");
  if (fw == NULL) goto error3;
  
  if (read_iges_line(fr, line) == -1) goto error5; /* end of file */
  if (line[72] == 'F') goto error6;       		   /* binary IGES */

  /* Skip S and G sections and find the 1st line of the 1st D entry,
   * then load 2nd line of the first D entry as well (2).
   * Analyze the pair of lines and if the enity is recognized, load it
   * into array entity[][] (3).
   * */
   
  while (line[72] != 'D') {
    err_code=read_iges_line(fr, line);
    if (err_code == -1) goto error5;
  }
  err_code=read_iges_line(fr, pairline);             /* (2) */
  if (err_code == -1) goto error5;
  parse_d_entry(line, pairline);                    /* (3) */

  /* Repeat the same on subsequent entries in D section */
  while (line[72] == 'D') {
    err_code=read_iges_line(fr, line);
    if (err_code == -1) goto error5;
    if (line[72] == 'D') {
      err_code=read_iges_line(fr, pairline);
      if (err_code == -1) goto error5;
      parse_d_entry(line, pairline);
      if (entity_sum >= MAXENTITY) goto error7;
    }
  }

  /* For every entity in array entity[] find PD (Parameter Data) line with the sequence nr
   * and process subsequent lines (find coordinates and write to text file).
   * */
  copy_field(line,73,7,pd_seq_s); pd_seq=atoi(pd_seq_s);


for (emark=0; emark < entity_sum; emark++) 
{

    while (pd_seq != entity[emark][PD_PTR]) {     /* find line in PD */
      read_iges_line(fr, line);
      copy_field(line,73,7,pd_seq_s); pd_seq=atoi(pd_seq_s);
    }
    if ( entity[emark][E_TYPE] == 126) 
    {
  
        if (o_write_layername == 1) 
            fprintf(fw, "%s\n", layer[entity[emark][ELAYER]]); 
        
        fprintf(fw, "\n");
        
        for (i=1; i <= entity[emark][PD_CNT]; i++) 
        {   
            
          // 126 means BSpline:
          read_pd126_writetxt(line, i, bSpline, splineId);
          read_iges_line(fr, line); 
          
          // do the necessary work to copy strings, etc.
          copy_field(line,73,7,pd_seq_s); pd_seq=atoi(pd_seq_s); 
        // runs until the end of file.
   	  
        }
        fprintf(fw, "\n");
        
        // increment and go to the next bspline:
        splineId++; 
        
    }

}

mesh->sizes.splines = splineId;

goto end;
 
error3: printf("error: File opening didn't succeeded!\n");
  goto end;
error4: printf("error: File closing didn't succeeded!\n");
  goto end;
error5: printf("error: IGES file is empty!\n");
  goto end;
error6: printf("error: IGES file is not of ASCII type!\n");
  goto end;
error7: printf("error: insufficient limit \"MAXENTITY\" for entity processing!\n");
  goto end;
end:
  return ;
}

