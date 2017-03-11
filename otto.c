/*      Copyright University Corporation for Atmospheric Research, 1995, 1998, 2010   */

/*
 *      otto:  calculate autocorrelation and power spectrum 
 */

	// use this to write out values for input, e.g., to R
#define CORR 0
	// use this to get old stock-forecasting function
#define STOCKS 1
	// use this for special iterative forecasting at high rate
#define ITERATING 0
#define EXTERN
#define XVIEW 1
#define GLADE 2
#define XML   3
#define GUI GLADE
#if(GUI == GLADE)
#include <glade/glade.h>
#endif
#if (GUI == XML)
#endif
#if (GUI == XVIEW)
                        // xview stuff removed now, expecting this to be dead
#endif
#include <gtk/gtk.h>
//
#include "lib/AlwaysNeeded.h"

#define DEBUG 0
#define BUFVARS 7
#define TWOPI 6.28318530717959
#define FIVETHIRDS 1.6666666666667
#define TWOTHIRDS 0.6666666666667
#define MAXSEG 8192
#define MAXSMOOTH 5000
#define XWINDOWS 0
#define PYPLT 1

#include "Xlib/ncg.struct"
struct NCG_parameters ncg;

struct Daphead1 dhd1;
struct Daphead1 dhd1out;      // added for Python interface
struct Dapdata  ddataout;
struct Dapmisc  miscout;


FILE           *fcstfile; 

 /* function declarations   */

#include "Xlib/Xanadu_function_declarations"
extern void    *GetMemory();
extern int      arc_();
extern int      ipcnt_();
extern int      rand();
extern void     clsgks_();
extern void     opngks_();
extern          plchhq_();
extern void     setusv_();
extern void     dashdb_();
extern void     ncar_();
extern void     title_();
extern void     fnote_();
extern void     idline_();
extern void     frstpt_(), vector_(), lablx_(), lably_(), frame_(), curve_();
extern void     head_();
double		Square(), Parzen(), Welch(), Hanning();
void            add_lscale();
void            quit();
void            plot();
void            master();
void            p_acv();
void            set_scales();
void            p_fft();
void            p_mem();
void            set_limits();
void            next();
void            last();
void            set_panel();
int             selected();
int             cselected();
int             wselected();
int             limits();
int             cspec();
void            metafile();
void            print_plot();
void            acv_process();
void            fft_process();
void            mem_process();
void            interact();
void            getdata();
void            plotdata();
void            plotfcst();
void            detrend();
void            acv();
void            ftran();
void            fftplt();
void            cosplt();
void            cophase();
void            cphase();
float           ftsec();
float           fsect();
int             itsec();
void            save_def();
void            pad();
int             kdate_();
float          *vector();



char            filename[MAX_FILE_NAME_LENGTH];	/* data file        */
char           *fname;
char            arcbuf[BUFVARS][NSZ] = {"WDC   ", "WSC   ", "THDG  ", "TASX  ", "WIC   ", "ATX   ", "TASX  "};
char            toplabel[80];

int             dbg = DEBUG;
int             iquit = 0;
int             date, titleline, today;
int             start_time, end_time, cycle_interval, accumulate;
int		addfft = 0, addmem = 0;
float           sample_frequency = 1.;
int             smooth_bins = 100;
float           duration, tstart;
int             data_size;	/* number of data elements */
int             buffer_size;	/* size of buffer containing data    */
int             fbuffer_size;	/* size of forecast buffer  */
float          *data_buffer;	/* storage for data, allocated later */
float          *fcst_buffer;	/* storage for forecast, allocated later */
float          *co_data_buffer;	/* storage for data, allocated later */
float          *cc_data_buffer;	/* storage for data, allocated later */
int             data_index;	/* pointer to data_buffer  */

int             label_length, ilbl;
int             ishowp = 0;
float           angd, cntr, size;	/* arguments for plchhq */
float           buf[BUFVARS + 2];
float           buff[2];
double         *r;		/* autocovariance array */
double          variance, flux;
int             max;
int             points, dpoints;
int             segment_length = 256;	/* must be power of 2 for fft */
int             poles = 50;
int             dash_pattern = 16711935, dash_pat2 = 252645135;
float           resolution = 0.01;
int		pday1 = 60, pday2 = 90;
char		labelF1[50], labelF2[50];
char            res_text[10];
char            spec_text[30];
float           vmin = 1.e30, vmax = -1.e30, tas_average = 0.;
float           cvmin = 1.e30, cvmax = -1.e30;
float           value;		/* variables stored, scr */
float		filtered_value = 0.;
float		filtered_5by3 = 0.;
float           cvalue;		/* variables stored, scr */
float		value_Last, value_Day1, value_Day2;
float		CurrentValue, CurrentCost, Profit;
double          mean = 0., trend = 0.;
double          co_mean = 0., co_trend = 0.;
double		ptrend;
int             window = 1, iwindow = 0;
int             show_errors = 0;
float           flow_spec, fhigh_spec, plow_spec, phigh_spec, fplow_spec, fphigh_spec,
		phlow_spec = -180., phhigh_spec = 180.,
                epslow_spec, epshigh_spec, cplow_spec, cphigh_spec;
float           flux_min_freq = 0., flux_max_freq = 100000.;
int             batfft = 1, batmem = 1;
			// pyplot controls generation of python routine for plot
			// print_plt controls batch-run send-to-printer
int             pyplot = 0, print_plt = 0;
char           *dapdata;
char            variable[NSZ]   = "WIC       ";	/* for spectrum */
char            covariable[NSZ] = "TASX      ";	/* for cospectrum */
char            ccvariable[NSZ] = "TASX      ";	/* for cost (forecast) */
char            wdvariable[NSZ] = "WDC       ";	/* for spectrum */
char            wsvariable[NSZ] = "WSC       ";	/* for spectrum */
char            hvariable[NSZ]  = "THDG      ";	/* for spectrum */
int		ishowfft, ishowacv, ishowmem, ishowflx, ishowcfft, ishowcacv, ishowcmem, ishowcflx;
int             spect_variable = 3;	/* index for button selection */
int             cospect_variable = 3;	/* index for button selection */
int             proj_degrees = 0;	/* angle to project wind along */
int             cproj_degrees = 0;	/* angle to project wind along */
int		delay = 0;		/* delay in ms		*/
int		idelay = 0;		/* delay in samples	*/
int             maxlags, maxlags_s = 600;
int             plotlags, plotlags_s = 600;
int             smoothpts, smooth_s = 150;
char            dummy = ' ';
char            par[8], dpar[8];
int             izero = 0, ione = 1, itwo = 2, ithree = 3, ifour = 4, ifive = 5,
                iseven = 7;
int             icolor = 2;
int             thousand = 1000, two_thousand = 2000, three_thousand = 3000;
float           mone = -1., zero = 0., one = 1.;
float           random_noise = 0.01, cosine_amplitude = 0.01, cosine_period = 5.5;
float		filter_tau = 0.;
float		random_5by3 = 0.;
float           co_random_noise = 0.01, co_cosine_amplitude = 0.01, co_cosine_period = 5.5;
float		co_cosine_phase = 0.0;
FILE           *temp_file;	/* scratch file */
double		LogRMS;
int		nxtf = 1;


int             restore();
void            show_dir();
void            setfile();
void            write_defaults();
int             fquit();


char            title_line[MAXLINE];
char            file_specs[MAXLINE];	/* plot specs            */
char           *drgv = "otto";

 /* include initial values for specs here */
int             end_def = 400000;

 /* plot specifications in spec file      */
struct specifications {
    char           *keyword;
    char            type;
    int            *variable_address;
    char           *window_item_name;
    int            *default_value;
}               spec[] = {

    "START", INT, &start_time, "START:", &izero,
    "END", INT, &end_time, "END:", &end_def,
    "CYCLE", INT, &cycle_interval, "Cycle interval:", &izero,
    "ADDTO", INT, &accumulate, "Accumulate", &accumulate,
    "ADDFFT", INT, &addfft, "Add fft", &addfft,
    "ADDMEM", INT, &addmem, "Add mem", &addmem,
    "SVAR", INT, &spect_variable, "VARIABLE", &spect_variable,
    "VAR", VAR, (int *) variable, "OTHER VARIABLE:", (int *) variable,
    "COVAR", VAR, (int *) covariable, "CO-VARIABLE:", (int *) variable,
    "WDVAR", VAR, (int *) wdvariable, "WIND DIR. VARIABLE:", (int *) wdvariable,
    "WSVAR", VAR, (int *) wsvariable, "WIND SPEED VARIABLE:", (int *) wsvariable,
    "HVAR", VAR, (int *) hvariable, "HEADING VARIABLE:", (int *) hvariable,
    "SHOWFFT", INT, &ishowfft, "showfft", &ithree,
    "SHOWACV", INT, &ishowacv, "showacv", &ithree,
    "SHOWMEM", INT, &ishowmem, "showmem", &ithree,
    "SHOWFLX", INT, &ishowflx, "showflx", &izero,
    "SHOWCFFT", INT, &ishowcfft, "showfft", &izero,
    "SHOWCACV", INT, &ishowcacv, "showacv", &izero,
    "SHOWCMEM", INT, &ishowcmem, "showmem", &izero,
    "SHOWCFLX", INT, &ishowcflx, "showflx", &izero,
    "PROJ", INT, &proj_degrees, "Projection for wind component (deg.):", &izero,
    "CPROJ", INT, &cproj_degrees, "Projection, wind component (deg.):", &izero,
    "DELAY", INT, &delay, "Delay (ms):", &delay,
    "MAXLAGS", INT, &maxlags_s, "Maximum lag for autocorrelation (s):", &maxlags_s,
    "MAXPLOT", INT, &plotlags_s, "Maximum lag to plot (s):", &plotlags_s,
    "SMOOTHS", INT, &smooth_s, "Smoothing time constant (s):", &smooth_s,
    "SEGL", INT, &segment_length, "Points per segment (power of 2):", &segment_length,
    "WINDOW", INT, &iwindow, "WINDOW for data segments:", &iwindow,
    "ERRORS", INT, &show_errors, "error bands?", &show_errors,
    "SMOOTHB", INT, &smooth_bins, "log smoothing intervals:", &smooth_bins,
    "POLES", INT, &poles, "Poles for Maximum-entropy method:", &poles,
    "RESN", FLT, (int *) &resolution, "Resolution for MEM:", (int *) &resolution,
    "PDAY1", INT, &pday1, "Forecast Days, #1:", &pday1,
    "PDAY2", INT, &pday2, "Forecast Days, #2:", &pday2,
    "LCOLOR", INT, &icolor, "  ", &icolor,
    "FLOW", FLT, (int *) &flow_spec, "Frequency lower limit:", (int *) &zero,
    "FHIGH", FLT, (int *) &fhigh_spec, "Frequency upper limit:", (int *) &zero,
    "PLOW", FLT, (int *) &plow_spec, "Spectrum lower limit:", (int *) &zero,
    "PHIGH", FLT, (int *) &phigh_spec, "Spectrum upper limit:", (int *) &zero,
    "WLOW", FLT, (int *) &fplow_spec, "Weighted spectrum lower limit:", (int *) &zero,
    "WHIGH", FLT, (int *) &fphigh_spec, "Weighted spectrum upper limit:", (int *) &zero,
    "ELOW", FLT, (int *) &epslow_spec, "Dissipation-rate spectrum lower limit:", (int *) &zero,
    "EHIGH", FLT, (int *) &epshigh_spec, "Dissipation-rate spectrum upper limit:", (int *) &zero,
    "CPLOW", FLT, (int *) &cplow_spec, "Cospectrum lower limit:", (int *) &zero,
    "CPHIGH", FLT, (int *) &cphigh_spec, "Cospectrum upper limit:", (int *) &zero,
    "PHLOW", FLT, (int *) &phlow_spec, "Phase lower limit:", (int *) &phlow_spec,
    "PHHIGH", FLT, (int *) &phhigh_spec, "Phase upper limit:", (int *) &phhigh_spec,
    "FXMINF", FLT, (int *) &flux_min_freq, "Minimum frequency for flux calc.", (int *) &flux_min_freq,
    "FXMAXF", FLT, (int *) &flux_max_freq, "Maximum frequency for flux calc.", (int *) &flux_max_freq,
    "BATFFT", INT, &batfft, "Batch FFT", &batfft,
    "BATMEM", INT, &batmem, "Batch MEM", &batmem,
    "PYPLOT", INT, &pyplot, "Python Plot", &izero,
    "PRINT",  INT, &print_plt, "Print plot, batch", &izero
};

int             number_specifications = sizeof(spec) / sizeof(struct specifications);

int             ipc = 0;
int             batch_mode = FALSE;

/*
 *	fileout routine: adapted from file.c, old RAF routine
 */
#define F_OK 0			/* including unistd.h gave error */
#define MAX_FLNM_LNGTH  80

char           *dapdata;	/* path for data file name */

fileout(int ioptn) 

     /*   DAP File Read/Write Subr.    MLG/RLR/CAW/MDD

	  Subroutine Arguments:
	  ioptn (input)  => Task to perform
	  4 => Write data (open data file, if necessary; leave open, exclusive)
	  5 => Close the data file, writing out the record first, if necessary)

	  ierr  (output) => Error return
	  0 => No error
	  1 => Illegal IOPTN
	  2 => Illegal record number
	  3 => Unexpected header EOF
	  <0 => FMP error number

	  lue   (input)  => Logical unit for reporting any error (0=none).


	  Assumption:  During program execution, the data file is intended
	  to be either read or written.  The data file is opened for shared
	  access by the first read request and opened for exclusive access
	  by the first write request.  It remains open in its shared/exclusive
	  mode for faster response to subsequent read/write requests.
	  */
{
  static int      fh;		/* header file descriptor */

  /* (data file is global)   */
  int             nbytes;
  int             result = -1;
  static char     ofilename[MAX_FLNM_LNGTH];
  char            headname[8], dataname[8];
  static int      fd_out, RecordCount;


  switch (ioptn) {

  case 4:	/****** Come here to write a data file record *****/
    if (!miscout.iopnfi) {
        system ("rm -f ./DXDATO");
		// count records output
        RecordCount = 0;
	strcpy (ofilename, "./DXDATO");
        printf(" output file name is %s\n", ofilename);
        if ((fd_out = open(ofilename, O_RDWR | O_CREAT | O_TRUNC, 0644)) == -1) {
	    perror("FILE.C(write data): error opening data file");
	    sleep(2);
	    exit(-1);
        }
        miscout.iopnfi = TRUE;
        printf(" open successful.\n");
    }
    if ((nbytes = write(fd_out, &ddataout, dhd1out.nwords * sizeof(float))) == -1) {
        perror("FILE.C(write data): write data");
        exit(-1);
    }
    RecordCount++;
    break;

  case 5:			/******** Come here to close the data file *******/
    close(fd_out);
    miscout.iopnfi = FALSE;
    break;

  }
  return (0);
}

void
write_defaults(dfile)
    char           *dfile;

{

    void read_all();
    int  i, j;

    specs_file = fopen(dfile, "w");
    if (specs_file == NULL) {
	return;
    }
    read_all();
    if ((int) strlen(title_line) > 1) {
	fprintf(specs_file, "#%s\n", title_line);
    }

    for (i = 0; i < number_specifications; i++) {
	switch (spec[i].type) {
	case INT:
	    if (!strcmp(spec[i].keyword, "START"))
		fprintf(specs_file, "START 0 \n");
	    else if (!strcmp(spec[i].keyword, "END"))
		fprintf(specs_file, "END 400000 \n");
	    else
		fprintf(specs_file, "%s %d\n", spec[i].keyword, *spec[i].variable_address);
	    break;
	case FLT:
	    fprintf(specs_file, "%s %f\n", spec[i].keyword, *((float *) spec[i].variable_address));
	    break;
	case VAR:
	    fprintf(specs_file, "%s %s\n", spec[i].keyword, (char *) spec[i].variable_address);
	    break;
	case FNM:
	case STR:
	case LIN:
	    fprintf(specs_file, "%s %s\n", spec[i].keyword, (char *) spec[i].variable_address);
	    break;
	case COL:
            fprintf(specs_file, "%s ", spec[i].keyword);
            for (j = 0; j < COLORED_ITEMS; j++)
                fprintf(specs_file, "%d ", (*(spec[i].variable_address + j) + 1));
            fprintf(specs_file, "\n");
            break;
        case VTR:
            if(tseg_max <= 0) break;
            fprintf(specs_file, "%s ", spec[i].keyword);
            for (j = 0; j < tseg_max; j++)
                fprintf(specs_file, "%8.2f", *((float *) (spec[i].variable_address + j)));
            fprintf(specs_file, "\n");
            break;
        case VTI:
            if(tseg_max <= 0) break;
            fprintf(specs_file, "%s ", spec[i].keyword);
            for (j = 0; j < tseg_max; j++)
                fprintf(specs_file, "%d", *(spec[i].variable_address + j));
            fprintf(specs_file, "\n");
            break;

	}
    }
		// special program-dependent code goes here
    fprintf(specs_file, "FIN\n");
    (void) fclose(specs_file);
}


int
restore(dfile)
    char           *dfile;

{

    char           *ch;

    char            line[MAXLINE + 1];
    char            par[NSZ];
    char            s1[20], s2[20];
    char            sysbuf[40];
    char	    tbuffer[10];
    int             n1, n2;
    int             i, j, k, m;
    float           low, high;
    static int      set_window = 0;
    short           x, y;
    char            pbuf[VNAME_SIZE + 6];


    title_line[0] = '\0';
    specs_file = fopen(dfile, "r");
    if (specs_file == NULL) {
	return (0);
    }
    /* set to start-up values */
    for (i = 0; i < number_specifications; i++) {
	switch (spec[i].type) {
	case INT:
	    *spec[i].variable_address = *spec[i].default_value;
	    break;
	case FLT:
	    *((float *) spec[i].variable_address) = *((float *) spec[i].default_value);
	    break;
	case VAR:
	    strncpy((char *) spec[i].variable_address, (char *) spec[i].default_value, NSZ);
	    PAD((char *) spec[i].variable_address);
	    printf(" restore set variable name to %s\n", (char *) spec[i].variable_address);
	    break;
	case FNM:
            strcpy((char *) spec[i].variable_address, (char *) getenv("DAPFILE"));
            break;
	case STR:
	case LIN:
	    strcpy((char *) spec[i].variable_address, (char *) spec[i].default_value);
	    break;
        case COL:
            colored_item = 0;
            for (j = 0; j < COLORED_ITEMS; j++)
                *(spec[i].variable_address + j) = *(spec[i].default_value + j) - 1;
            break;
        case VTR:
            tseg = 0;
            tseg_max = 0;
            for (j = 0; j < NTSEGS; j++)
                *((float *)spec[i].variable_address + j) = *((float *)spec[i].default_value + j) - 1;
            break;
        case VTI:
            tseg = 0;
            tseg_max = 0;
            for (j = 0; j < NTSEGS; j++)
                *(spec[i].variable_address + j) = *(spec[i].default_value + j) - 1;
            break;

	default:
	    break;
	}
    }
//    printf ("before restore, icolor=%d\n", icolor);
    /* read spec file and set values */
    while (fgets(line, MAXLINE, specs_file) != NULL) {
	if (!strncmp(line, "FIN", 3))
	    break;
	sscanf(line, "%s", s1);
	i = -1;
	for (j = 0; j < number_specifications; j++) {
	    if (!strcmp(s1, spec[j].keyword)) { i = j; }
	}
	if (i != -1) {
	    switch (spec[i].type) {
	    case INT:
		sscanf(line, "%s %d", s2, spec[i].variable_address);
  		printf(" setting integer spec %d to %d\n", i, *spec[i].variable_address);
		break;
	    case FLT:
		sscanf(line, "%s %f", s2, (float *) spec[i].variable_address);
		printf(" setting float spec %d to %f\n", i, *(float *) spec[i].variable_address);
		break;
	    case VAR:
		sscanf(line, "%s %s", s2, (char *) spec[i].variable_address);
	        PAD((char *) spec[i].variable_address);
	        printf(" restore set variable name to %s\n", (char *) spec[i].variable_address);
		break;
            case FNM:
                sscanf(line, "%s %5s", s2, (char *) spec[i].variable_address);
                for(j=0;j<5;j++) {
                        *((char *)spec[i].variable_address+j) =
                            toupper(*((char *)spec[i].variable_address+j));
                        if(!isalnum(*((char *)spec[i].variable_address+j)) &&
                                        *((char *)spec[i].variable_address+j) !=
 BLANK) {
                            for(k=j;k<6;k++)
                                *((char *)spec[i].variable_address+k) = BLANK;
                            break;
                        }
                }

	    case STR:
		sscanf(line, "%s %s", s2, (char *) spec[i].variable_address);
		break;
            case LIN:
                strncpy((char *) spec[i].variable_address, (line+8), 80);
                break;
            case COL:
                ch = strtok(line, WHITE);
                j = 0;
                while (ch != NULL && j < COLORED_ITEMS) {
                    if ((ch = strtok(NULL, WHITE)) == NULL)
                        break;
                    if (sscanf(ch, "%d", (spec[i].variable_address + j)) == EOF)
                        break;
                    (*(spec[i].variable_address + j++))--;
                    k = j - 1;
                }
                colored_item = 0;
                break;
            case VTR:
                ch = strtok(line, WHITE);
                j = 0;
                while (ch != NULL && j < NTSEGS) {
                    if ((ch = strtok(NULL, WHITE)) == NULL)
                        break;
                    if (sscanf(ch, "%f", ((float *)spec[i].variable_address + j)) == EOF)
                        break;
                    (*((float *)spec[i].variable_address + j++))--;
                    k = j - 1;
                }
                tseg = 0;
                tseg_max = j;
                break;
            case VTI:
                ch = strtok(line, WHITE);
                j = 0;
                while (ch != NULL && j < NTSEGS) {
                    if ((ch = strtok(NULL, WHITE)) == NULL)
                        break;
                    if (sscanf(ch, "%d", (spec[i].variable_address + j)) == EOF)
                        break;
                    (*(spec[i].variable_address + j++))--;
                    k = j - 1;
                }
                tseg = 0;
                tseg_max = j;
                break;

	    }
	} else if (line[0] == '#') {	//provide for comments
	    if (title_line[0] == '\0') {
	        title_line[0] = '#';
	        strncpy(&(title_line[1]), &line[1], MAXLINE-1);	/* provide for comments */
	        if ((j = strlen(title_line)) > 0)
		    title_line[j - 1] = '\0';
	    }
	} 
    }				/* end of "while" loop to read from otto.def		 */
//    printf (" after restore, icolor=%d\n", icolor);
    (void) fclose(specs_file);
    set_window = 1;
    return (1);
}

#include "OttoGTKcallbacks.c"

int
main(int argc, char **argv)

{
    /* variable declarations   */
    int             i, j, k, l, m, n, nbytes;
    int             stime = 0, nvars = 1;
    int             brgc;
    char           *brgv[30];
    char            line[MAXLINE + 1];
    char            sample_f_buf[12];
    static char     gks[50];

    /* start of executable code */

    if (argc > 1)
	batch_mode = TRUE;
    fname = (char *) getenv("XANFILE");
    if ((fname == 0) || (*fname == ' ') || (*fname == '\0')) {
	printf(" DAP file name? ");
	scanf("%5s", filename);
    } else
	strncpy(filename, fname, MAX_FILE_NAME_LENGTH);

    start_time = 0;
    end_time = 400000;

    sprintf(gks, "NCARG_GKS_OUTPUT=%s", (char *) getenv("NCARG_GKS_OUTPUT"));
    strcat(gks, ".otto");
    putenv(gks);

    gtk_init(&argc, &argv);
    gtk_rc_parse ("./.XanaduGTK.rc");
    if (batch_mode) {
    } else {
#if (GUI == GLADE)
    char *glade_file;
    glade_file = (char *) malloc(150 * sizeof(char));
    strcpy(glade_file, (char *) getenv("XANADU"));
    strcat(glade_file, "/gui/Otto.glade");
    xml = glade_xml_new(glade_file, NULL, NULL);
    free(glade_file);
    set_GTK_callback_references();
#endif
#if (GUI == XML)
    char *xml_file;
    xml_file = (char *) malloc(150 * sizeof(char));
    strcpy(xml_file, (char *) getenv("XANADU"));
    strcat(xml_file, "/gui/Otto.xml");
    builder = gtk_builder_new ();
    gtk_builder_add_from_file (builder, xml_file, NULL);
    free(xml_file);
//  gtk_builder_connect_signals (builder, NULL);
    set_GTK_callback_references();
    g_object_unref (G_OBJECT (builder));
#endif
    }
//

    /* first call to fill in dhd1.ideltt */
    if ((j = arc_(filename, &nvars, buf, &stime, arcbuf)) < 0) {
	fprintf(stderr, "otto: unable to open data file\n");
	printf(" otto: unable to open data file\n");
	return;
    }
#if(0)
    printf("dhd1.ideltt=%d, file_type=%c, proj=%d,date=%d,words=%d,segs=%d\n",
	dhd1.ideltt, dhd1.file_type, dhd1.iproj, dhd1.idated, dhd1.nwords,
	dhd1.ntmseg);
#endif
    sample_frequency = 1000. / dhd1.ideltt;


    /*
     * set time limits, variable, etc, then call "process" for calculations
     */

    strcpy(file_specs, "otto.def"); 
    (void) restore(file_specs);
    if (!strncmp("LONW", variable, 4))
	spect_variable = 0;
    else if (!strncmp("LATW", variable, 4))
	spect_variable = 1;
    else if (!strncmp("WIND", variable, 4))
	spect_variable = 2;
    strncpy(arcbuf[0], wdvariable, NSZ);
    strncpy(arcbuf[1], wsvariable, NSZ);
    strncpy(arcbuf[2], hvariable, NSZ);
    window = iwindow + 1;
    if (batch_mode) {
	system("rm -f $NCARG_GKS_OUTPUT");
	maxlags = maxlags_s * sample_frequency;
	plotlags = plotlags_s * sample_frequency;
	smoothpts = smooth_s * sample_frequency;
	if (batfft) {
	    fft_process();
    	    clsgks_();
	    if (print_plt) {system(SEND_COMMAND);}
	    system("mv -f ./$NCARG_GKS_OUTPUT $NCARG_GKS_OUTPUT.fft");
	}
	if (batmem) {
#ifndef NR
	    fprintf(stderr, " MEM analysis not supported, this version\n");
	    fprintf(stderr, " (MEM spectral analysis needs Numerical Recipes routines)\n");
#else
//            printf (" before mem_process, icolor=%d\n", icolor);
	    mem_process(addmem, icolor);
    	    clsgks_();
  	    if (print_plt) {system(SEND_COMMAND);}
	    system("mv -f ./$NCARG_GKS_OUTPUT $NCARG_GKS_OUTPUT.mem");
#endif
	}
	exit(0);
    }
    (void) ResetPanel();
    brgv[0] = "variance spectra";
    brgc = argc;
    if (brgc > 28)
	brgc = 28;
    for (i = 0; i < brgc; i++)
	brgv[i] = argv[i];
    brgv[brgc] = "-geometry";
    brgv[brgc + 1] = "1150x220+0-0";
    brgc += 2;
    brgv[brgc] = NULL;
    maxlags = maxlags_s * sample_frequency;
    plotlags = plotlags_s * sample_frequency;
    smoothpts = smooth_s * sample_frequency;
    (void) sprintf(res_text, "%8.3f", resolution);

    gtk_main();

    exit(0);
}


void
acv_process(accumulate)
    int             accumulate;

{
    int		    nv = 4, stime;
    opngks_();
				/* determine variables needed, load and call
				 * arc to set ideltt appropriate for those
				 * variables. */
    if (!strncmp(variable, "LONW", 4) || !strncmp(variable, "Long",4)) {
    } else if (!strncmp(variable, "LATW", 4) || !strncmp(variable, "Late",4)) {
    } else if (!strncmp(variable, "WIND", 4) || !strncmp(variable, "Proj",4)) {
    } else if (!strncmp(variable, "SIM ", 4) || !strncmp(variable, "Simu",4)) {
    } else {
	strncpy(arcbuf[4], variable, NSZ);
	nv = 5;
    }
    strcpy(arcbuf[0], wdvariable);
    strcpy(arcbuf[1], wsvariable);
    strcpy(arcbuf[2], hvariable);
    strcpy(arcbuf[3], "TASX  ");
    stime = 0;
				/* this sets ideltt ... */
    (void) arc_(filename, &nv, buf, &stime, arcbuf);
    sample_frequency = 1000. / dhd1.ideltt;

    duration = ftsec((float) end_time) - ftsec((float) start_time);
    tstart = ftsec((float) start_time);
    idelay = delay / dhd1.ideltt;
    if (idelay < 0)
	idelay = 0;

    getdata(0);
    if (ishowacv & 1)
	plotdata(data_buffer, mean, trend, variable, proj_degrees, &vmin, &vmax);
    detrend(data_buffer, mean, trend, 1);
    acv(accumulate);
    if (accumulate)
	return;
    ftran();
    return;
}
void
fft_process()
{
    void            spctrm();
    double          (*win)();


    int             i, m, k, ovrlap = 1;
    int		    mi, ki;
    int		    nbad = 0;
    int		    nv = 4;
    int		    stime;
    float	    xtmp;
    double          variance_factor, pad_weight;
    double         *fps, *fpsc, *p, *q, *dw;

    opngks_();
				/* determine variables needed, load and call
				 * arc to set ideltt appropriate for those
				 * variables. */
    printf(" now in fft_process, ishowfft=%d, ishowcfft=%d, variable=(%s), covariable=(%s)\n", ishowfft, ishowcfft, variable, covariable);
    if (!strncmp(variable, "LONW", 4)) {
    } else if (!strncmp(variable, "LATW", 4)) {
    } else if (!strncmp(variable, "WIND", 4)) {
    } else if (!strncmp(variable, "SIM ", 4)) {
    } else {
	strncpy(arcbuf[4], variable, NSZ);
	nv = 5;
    }
    strcpy(arcbuf[0], wdvariable);
    strcpy(arcbuf[1], wsvariable);
    strcpy(arcbuf[2], hvariable);
    strcpy(arcbuf[3], "TASX      ");
    printf(" in fft_process, variables (WD,WS,H) are (%s) (%s) (%s)\n",
	   wdvariable, wsvariable, hvariable);
    if (ishowcfft) {
	if (!strncmp(covariable, "LONW", 4)) {
	} else if (!strncmp(covariable, "LATW", 4)) {
	} else if (!strncmp(covariable, "WIND", 4)) {
	} else if (!strncmp(covariable, "SIM ", 4)) {
	} else if (!strncmp(covariable, "NONE", 4)) {
	} else {
	    strncpy(arcbuf[5], covariable, NSZ);
	    nv = 6;
	}
    }
    stime = 0;
				/* this sets ideltt ... */
    (void) arc_(filename, &nv, buf, &stime, arcbuf);
    sample_frequency = 1000. / dhd1.ideltt;

    duration = ftsec((float) end_time) - ftsec((float) start_time);
    tstart = ftsec((float) start_time);
    idelay = delay / dhd1.ideltt;
    if (idelay < 0)
	idelay = 0;

    getdata(1);
    if (ishowfft & 1)
	plotdata(data_buffer, mean, trend, variable, proj_degrees, &vmin, &vmax);
    if (ishowcfft & 1)
	plotdata(co_data_buffer, co_mean, co_trend, covariable, cproj_degrees,
		 &cvmin, &cvmax);
			// 0 because missing_v set == 0 after variance calc.
    detrend(data_buffer, mean, trend, 0);
				/* calculate total variance, time series */
    variance = 0.;
    for (i = 0; i < buffer_size; i++) {
	if (*(data_buffer + i) != MISSING_DATA) 
	    variance += *(data_buffer + i) * *(data_buffer + i);
	else {
	    nbad++;
	    *(data_buffer + i) = 0.;
	}
    }
    variance /= (double) (data_size - nbad);
    printf(" variance from time series (after mean and trend removed) = %f\n",
	(float) variance);
    m = segment_length / 2;
		/* was k = points / segment_length - 1; */
    k = (int) (xtmp = ((float)data_size / (float) segment_length) - 0.5);
    if ((float) k < xtmp) k++;
   // fps = (double *) GetMemory((size_t) (m+1) * sizeof(double));
    printf(" attempt #1 to allocate %d double words\n", m+1);
    fps = (double *) GetMemory((size_t) (m+1) * sizeof(double));
    if (k > 0) {
#if(0)
	printf(" m,k,points=%d %d %d\n");
	for (ki = 0; ki < k; ki++)
	    for(mi = 0; mi < 2*m; mi++) {
		if (mi%10 == 0) printf(" %d %d ", ki, mi);
		printf("%f ", *(data_buffer+mi+ki*segment_length));
		if (mi%10 == 9) printf("\n");
	    }
#endif
	window = iwindow + 1;
#ifndef NR		/* define in Makefile, via -DNR; without this, uses
			 * Chris Webster's versions of fft.c and spctrm */
	switch (window) {
	case 1:		/* Parzen	*/
	  win = Parzen;
	  break;
	case 2:		/* Square	*/
	  win = Square;
	  break;
	case 3:		/* Welch	*/
	  win = Welch;
	  break;
	case 4:		/* Hanning	*/
	  win = Hanning;
	  break;
	}
	Spectrum(data_buffer, fps, 2 * k, m, win);
#else
	printf(" call to spctrm with m=%d, window=%d\n", m, window);
	spctrm(data_buffer, fps, m, k,
	       ovrlap, window);
#endif
    } else {
	printf(" too few data points for segmenting with selected segment length\n");
	fprintf(stderr, " too few data points for selected segment length\n");
	return;
    }
				/* correct for zero padding: */
				/* (but don't use points; already padded...*/
    pad_weight = (double)((2*k+1)*m) / ((double) data_size);
printf(" pad_weight=%f, k=%d, m=%d, data_size=%d, buffer_size=%d\n",
         (float)pad_weight, k, m, data_size, buffer_size);

				/* need to multiply fps by this factor */
    variance_factor = 9. * k / 11.;
    if (variance_factor < 1.)
	variance_factor = 1.;
				/* plot the variance spectrum */
    fftplt(fps, variance_factor, pad_weight, spect_variable, variable, proj_degrees, ishowfft);
    if (ishowcfft) {
	detrend(co_data_buffer, co_mean, co_trend, 1);
	if (k > 0) {
    printf(" attempt #2 to allocate %d double words\n", m+1);
	    fpsc = (double *) GetMemory((unsigned) (m+1) * sizeof(double));
#ifndef NR
	Spectrum(co_data_buffer, fpsc, 2 * k, m, win);
#else
	    spctrm(co_data_buffer, fpsc, m, k,
		   ovrlap, window);
#endif
	} else {
	    printf(" too few data points for selected segment length\n");
	    fprintf(stderr, " too few data points for selected segment length\n");
	    return;
	}
	fftplt(fpsc, variance_factor, pad_weight, cospect_variable, covariable, cproj_degrees, ishowcfft);
	if (ishowcfft &  48) {
    printf(" attempt #3  to allocate %d double words\n", m+1);
	    p = (double *) GetMemory((size_t) (m+1) * sizeof(double));
    printf(" attempt #4 to allocate %d double words\n", m+1);
	    q = (double *) GetMemory((size_t) (m+1) * sizeof(double));
#ifndef NR
	    CoSpectrum(data_buffer, co_data_buffer, p, q, 2 * k, m, win);
#else
//	    cospec(data_buffer, co_data_buffer, p, q, m, k, ovrlap, window);
//	temporary: recalculate fps and fpsc in cospec (testing):
  	    cospec(data_buffer, co_data_buffer, p, q, fps, fpsc, m, k, ovrlap, window);
#endif
#if(1)
            printf(" start of cospec results\n");
	    for (i = 0; i < m; i++) {
		printf(" %d %e %e %e %e %lf\n",
		         i, fps[i]*pad_weight, fpsc[i]*pad_weight, p[i]*pad_weight, q[i]*pad_weight, (p[i]*p[i]+q[i]*q[i])/(fps[i]*fpsc[i]));
	    }
            printf(" end of cospec results\n");
            printf(" variance factor is %f\n", variance_factor);
#endif
	    if (ishowcfft & 16) {
		cosplt(p, q, variance_factor, spect_variable, variable,
		   cospect_variable, covariable, proj_degrees,
		   cproj_degrees);
	    }
	    if (ishowcfft & 32) {
		cophase(p, q, fps, fpsc, variance_factor, spect_variable,
			variable, cospect_variable, covariable, proj_degrees,
			cproj_degrees);
	    }
	}
    }
    free(fps);
    free(data_buffer);
    if (ishowcfft) {
	free(fpsc);
	if (ishowcfft & 48) {
	    free(p);
	    free(q);
	}
	free(co_data_buffer);
    }
    return;
}

void
mem_process(addmem, jcolor)
    int	addmem;
    int jcolor;
{
    float           evlmem();
    void            reimem();
    int             ipc, ipcold;
    float           pm, pmc, pmcv[4];
    float          *cof, *cofc, *cofcv;
    int		    status_mcar;
    float           frequency;
    float          *ps, *psc, *nu;
    float          *cospectrum, *quadrature;
    float           upl, botl;
    float           flow = .001;
    float           fhigh = 10.;
    float	    real, imag;
    float           a, a1, a2;
    double          b, c1, c2, c3, c4;
    float           re1, re2, im1, im2;
    float           df, logf;
    float           fdt;
    float           last_f;
    int             max2;
    int             label_length;
    char            label[46];
    float           fmin, fmax;
    float           variance_20km;	/* special KOFSE  */
    double          variance_factor;
    int             m_eps;
    float           eps, eps_bar, eps2_bar;
    float           ae;
    float	    t;
    double          sigma;
    double	    theta;
    int             isgn, ifrst;

    int             i, j, k, l, m, n;
    int		    data_index;
    int             iset;
    int		    nv = 4, stime;
    double	    localRMS;

    opngks_();

				/* determine variables needed, load and call
				 * arc to set ideltt appropriate for those
				 * variables. */
    printf (" in mem_process, jcolor=%d\n", jcolor);
    if (!strncmp(variable, "LONW", 4)) {
    } else if (!strncmp(variable, "LATW", 4)) {
    } else if (!strncmp(variable, "WIND", 4)) {
    } else if (!strncmp(variable, "SIM ", 4)) {
    } else {
	strncpy(arcbuf[4], variable, NSZ);
	nv = 5;
    }
    printf ("nv=%d, arcbuf[4] is %s\n", nv, arcbuf[4]);
    strcpy(arcbuf[0], wdvariable);
    strcpy(arcbuf[1], wsvariable);
    strcpy(arcbuf[2], hvariable);
    strcpy(arcbuf[3], "TASX  ");
    if (ishowcmem) {
	if (!strncmp(covariable, "LONW", 4)) {
	} else if (!strncmp(covariable, "LATW", 4)) {
	} else if (!strncmp(covariable, "WIND", 4)) {
	} else if (!strncmp(covariable, "SIM ", 4)) {
	} else {
	    strncpy(arcbuf[5], covariable, NSZ);
	    nv = 6;
	}
    }
#if(STOCKS)
    if ((ishowmem & 16)) {	/* forecast plot */
	strncpy(arcbuf[5], variable, NSZ);
	arcbuf[5][0]='N';
	strncpy(arcbuf[6], variable, NSZ);
	arcbuf[6][0]='C';
	nv = 7;
	ishowp = 1;
    }
#endif
    stime = 0;
				/* this sets ideltt ... */
    (void) arc_(filename, &nv, buf, &stime, arcbuf);
    sample_frequency = 1000. / dhd1.ideltt;

    tstart = ftsec((float) start_time);
    duration = ftsec((float) end_time) - tstart;
    idelay = delay / dhd1.ideltt;
    if (idelay < 0)
	idelay = 0;

    getdata(0);
    printf ("data [0], [1], ... [5] %f %f %f %f %f %f\n",
	*(data_buffer),
	*(data_buffer+1),
	*(data_buffer+2),
	*(data_buffer+3),
	*(data_buffer+4),
	*(data_buffer+5));
    printf(" returned from getdata; ready for plotdata\n");
    if (ishowmem & 1) {
	plotdata(data_buffer, mean, trend, variable, proj_degrees, &vmin, &vmax);
    }
    printf(" after plotdata for first variable\n");
    if (ishowcmem & 1) {
	plotdata(co_data_buffer, co_mean, co_trend, covariable, cproj_degrees,
		 &cvmin, &cvmax);
    }
    printf(" after plotdata for both variables\n");
    cof = vector(1, poles+2);	/* allocate space for coefficients */
    max2 = 1. / resolution;
    ps = vector(1, max2+2);	/* and for variance spectrum */
    nu = vector(1, max2+2);
    detrend(data_buffer, mean, trend, 1);
    if (ishowcmem) {
	detrend(co_data_buffer, co_mean, co_trend, 1);
	cofc = vector(1, poles);
	psc = vector(1, max2);
    }
    sigma = ((duration * 1000. / dhd1.ideltt) - poles - max2 - 2);
    if (sigma > 0.)
	sigma = sqrt(sigma);
    else
	sigma = 1.1;
    for (iset = 0; iset < 2; iset++) {
	flow = 0.001;
	fhigh = 10.;
				/* the NR switch is for compilation/load
				 * purposes; routine should never be
				 * called if NR not defined...           */
    if (iset == 0) {
        fcst_buffer = NULL;
        printf(" opening forecast.text for writing\n");
        fcstfile = fopen("forecast.text", "w");
    }
#if(ITERATING)
    int iterator = 1000 / dhd1.ideltt;
    do {
#endif
#ifdef NR
	if (iset == 0) {
	    memcof(data_buffer-1, data_size, poles, &pm, cof);
	    printf(" have memcof for first variable; pm=%f, data_size=%d poles=%d\n", pm,data_size,poles);
	} else if (ishowcmem <= 1) {
	    continue;
	} else {
	    memcof(co_data_buffer-1, data_size, poles, &pmc, cofc);
	}
#endif

	printf(" have memcof for both variables\n");
	if (flow_spec > 0.)
	    flow = flow_spec;
	if (fhigh_spec > 0.)
	    fhigh = fhigh_spec;
	fmin = 1. / points;
//	fmin = 0.00005;
	fmax = 0.5;
	fmin = log((double) fmin);
	fmax = log((double) fmax);
	df = (fmax - fmin) / max2;
	fdt = 0.0;
	variance = 0.0;
	variance_20km = 0.;	/* special, KOFSE analysis */
	for (i = 0; i < max2; i++) {
	    last_f = fdt;
	    fdt = fmin + df * (i+0.5);
	    fdt = exp((double) fdt);
	    if (last_f == 0.) {last_f = fdt;} // skip first point

	    /*
	     * change convention to match PSD from 0 to fmax (factor of 2
	     * here)
	     */
#ifdef NR
	    if (iset == 0) {
		ps[i] = 2. * evlmem(fdt, cof, poles, pm) / sample_frequency;
	    } else {
		psc[i] = 2. * evlmem(fdt, cofc, poles, pmc) / sample_frequency;
	    }
#endif
	    nu[i] = fdt * sample_frequency;
	    if (iset == 0)
		variance += ps[i] * (fdt - last_f);
	    else
		variance += psc[i] * (fdt - last_f);
#if(0)
	    if (nu[i] < 0.5) {
	        printf(" variance contribution from freq %f is %f, ps=%f\n", nu[i],ps[i]*(fdt-last_f), ps[i]);
	    }
#endif
	    if (nu[i] > tas_average / 20000.) {
		if (iset == 0)
		    variance_20km += ps[i] * (fdt - last_f);
		else
		    variance_20km += psc[i] * (fdt - last_f);
	    }
	    ipc = i * 40 / max2 + 60;
	    if (ipc != ipcold) {
		ipcold = ipc;
	    }
	}
	/* now omit this; see above */
	/* variance *= 2.; */
	variance *= sample_frequency;
	variance_20km *= sample_frequency;

	if (sample_frequency > 10.) {
	    fhigh *= 10.;
	    botl /= 100.;
	}
	if (((iset == 0) && (ishowmem & 2)) || ((iset == 1) && (ishowcmem & 2))) {	/* P(f) */
	    i = log10((double) variance) + 3.;
	    upl = pow((double) 10., (double) i);
	    botl = pow((double) 10., (double) (i - 6));
	    if (upl > 1.e20 || botl < 1.e-20) {
		printf(" plot outside expected limits, low=%f, high=%f\n",
		       botl, upl);
		return;
	    }
	    if (plow_spec > 0.)
		botl = plow_spec;
	    if (phigh_spec > 0.)
		upl = phigh_spec;
	    ncar_(&flow, &fhigh, &botl, &upl, &ione, &ione, &ione, &ione, &ifour, &ione);
	    lablx_(":PGL:T:RU: :L:(:U:H:L:Z):U:", 27);
	    ilbl = 0;
	    if (iset == 0) {
	        if (!strncmp(variable, "WI", 2)) ilbl = 1;
	        if (!strncmp(variable, "LONW", 4)) ilbl = 1;
	        if (!strncmp(variable, "LATW", 4)) ilbl = 1;
	        if (!strncmp(variable, "LATW", 4)) ilbl = 1;
	    } else {
	        if (!strncmp(covariable, "WI", 2)) ilbl = 1;
	        if (!strncmp(covariable, "LONW", 4)) ilbl = 1;
	        if (!strncmp(covariable, "LATW", 4)) ilbl = 1;
	        if (!strncmp(covariable, "LATW", 4)) ilbl = 1;
	    }
	    if (ilbl)
	        lably_("P(:GL:T:RU:) [m:S1:2 s:S2:-1  ]", 31);
	    else
	        lably_("P(:GL:T:RU:)", 12);
	    add_lscale(&flow, &fhigh, &botl, &upl);
	    gcolor_(&icolor);

	    /* now average together values within smooth_bins */
	    for (isgn = -1; isgn <= 1; isgn++) {
		if (isgn && !show_errors)
		    continue;
		if (isgn) {
		    gcolor_(&ithree);
		    setusv_("LW", &thousand, 2);
		} else {
		    gcolor_(&icolor);
		    setusv_("LW", &two_thousand, 2);
		}
		df = (log10((double) fhigh) - log10((double) flow)) / smooth_bins;
		logf = log10((double) nu[1]);
//		frstpt_(nu, &a);
		ifrst = 1;
		i = 0;
		a = 0.;
		frequency = 0.;
		m = 0;
		while (i < max2-1) {	// change 110803
		    i++;
		    while (log10((double) nu[i]) - logf > df) {
			if (m != 0) {
			    a /= m;
			    a *= (1. + isgn / (sigma * sqrt((double) m)));
			    frequency /= m;
			    if (frequency < flow) {ifrst = 1.;}
			    if (a < botl)
				a = botl;
			    if (ifrst) {
				ifrst = 0;
				frstpt_(&frequency, &a);
			    } else
				vector_(&frequency, &a);
			    a = 0.;
			    frequency = 0.;
			    m = 0;
			}
			logf += df;
		    }
		    frequency += nu[i];
		    if (iset == 0)
			a += ps[i];
		    else
			a += psc[i];
		    m++;
		    if (smooth_bins == MAXSMOOTH) {
			a *= (1. + isgn / sigma);
			if (ifrst) {
			    ifrst = 0;
			    frstpt_(&frequency, &a);
			} else
			    vector_(&frequency, &a);
			frequency = a = 0.;
			m = 0;
		    }
		}
		if (m != 0) {
		    a /= m;
		    a *= (1. + isgn / (sigma * sqrt((double) m)));
		    frequency /= m;
		    if (a < botl)
			a = botl;
		    vector_(&frequency, &a);
		}
	    }			/* end of isgn loop */
	    gcolor_(&ione);
	    a = fhigh * 0.9;
	    upl *= 0.5;
	    cntr = 1.;
	    size = 12.;
	    angd = 0.;
	    (void) sprintf(label, "VARIANCE=%15.7g", variance);
	    label_length = 24;
	    /* plchhq_(&a, &upl, label, &size, &angd, &cntr, label_length); */
	    upl *= 2.;
	    head_(&date, &start_time, &end_time, &dummy, 1);
	    titleline = 1;
	    if (iset == 0)
		fnote_(&titleline, variable, strlen(variable));
	    else
		fnote_(&titleline, covariable, strlen(covariable));
	    idline_(&titleline, label, 24);
	    titleline = 2;
	    fnote_(&titleline, "MEM", 3);
	    (void) sprintf(label, "poles=%5d, resln=%8.6f, smooth bins=%5d",
			   poles, resolution, smooth_bins);
	    idline_(&titleline, label, 46);
	    if (iset == 0) {
		if (!strncmp(variable, "WIND", 4)) {
		    titleline = 3;
		    sprintf(label, "wind in %d degree direction", proj_degrees);
		    idline_(&titleline, label, strlen(label));
		}
	    } else {
		if (!strncmp(covariable, "WIND", 4)) {
		    titleline = 3;
		    sprintf(label, "wind in %d degree direction", cproj_degrees);
		    idline_(&titleline, label, strlen(label));
		}
	    }
	    gcolor_(&iseven);
	    gsclip_(&ione);
	    fhigh /= 10.;
	    a = botl * pow((fhigh / flow), 1.666667);
	    line_(&flow, &a, &fhigh, &botl);
	    fhigh *= 10.;
	    gcolor_(&ione);

	    /*
	     * this strange call forces color reset, needed if plot frames
	     * are overlaid by 'med'
	     */
	    line_(&thousand, &thousand, &thousand, &thousand);
	    frame_();
	}
	printf(" finished with P(f) section\n");
	if (((ishowmem & 4) && (iset == 0)) || ((ishowcmem & 4) && (iset == 1))) {/* f P(f) */
	    upl = 1.e-10;
	    botl = 1.e10;
	    eps_bar = 0.;
	    eps2_bar = 0.;
	    m_eps = 0;
	    if (iset == 0) {
		if (spect_variable == 0 || !strncmp (variable, "UX", 2)) {
		    spect_variable = 0;
		    ae = 1.5;
		} else 
		    ae = 2.0;
	    } else {
		if (cospect_variable == 0 || !strncmp (covariable, "UX", 2)) {
		    spect_variable = 0;
		    ae = 1.5;
		} else
		    ae = 2.0;
	    }
	    ae = pow(ae, (double) 1.5) * TWOPI / tas_average;
	    for (i = 1; i < max2; i++) {	// change 110805
		/* get eddy diss. rate */
// corrected Oct 2013: take average first, then 3/2 power:
//		if (nu[i] >= nu[max2-1] / 10.) {
  		if ((nu[i] > 0.1) && (nu[i] <= 8.)) {
		    if (iset == 0) {
//			a = pow((double) ps[i], (double) (1.5)) * ae
//			    * pow((double) nu[i], 2.5);
			a = ps[i] * pow ((double) nu[i], FIVETHIRDS);
		    } else {
//			a = pow((double) psc[i], (double) (1.5)) * ae
//			    * pow((double) nu[i], 2.5);
			a = psc[i] * pow ((double) nu[i], FIVETHIRDS);
		    }
		    eps_bar += a;
		    eps2_bar += a * a;
		    m_eps++;
		}
	    }
	    if (m_eps > 0) {
		eps_bar /= m_eps;
		eps2_bar /= m_eps;
		eps2_bar = sqrt((eps2_bar - eps_bar * eps_bar));
		eps_bar = ae * pow ((double) eps_bar, 1.5); // corr Oct 2013
		eps2_bar = ae * pow ((double) eps2_bar, 1.5);  // corr Oct 2013
		printf(" MEM average eddy dissipation rate=%.5e +/- %.5e\n",
		       eps_bar, eps2_bar);
	    }
		// write a separate python data file, same format as for FFT:
            float t = 0.;
            dhd1out.nwords = 7;
	    for (i = 1; i < max2; i++) {
                t += 0.001 * dhd1.ideltt;
                ddataout.ihr = (int) (t+0.001) / 3600;
                ddataout.imin = (int) ((t+0.001)-3600.*ddataout.ihr) / 60;
                ddataout.isec = ((int) (t+0.001)) % 60;
                ddataout.imsec = 0;
                ddataout.values[0] = t;
                ddataout.values[1] = nu[i];
                ddataout.values[2] = ps[i] * nu[i];
                (void) fileout (4);
            }
            (void) fileout (5);	// close the file
            char syscmd[50];	// and rename it to avoid overlap
            if (addmem) {
		nxtf = 1;
		jcolor += 1;
        	char *rdbuf;
        	rdbuf = (char *) malloc(121 * sizeof(char));
		sprintf (rdbuf, "DXDADD%d", nxtf);
		while (temp_file = fopen (rdbuf, "r")) {
		    fclose (temp_file);
		    nxtf++;
		    jcolor += 1;
		    sprintf (rdbuf, "DXDADD%d", nxtf);
		    printf ("ADD file exists; advancing nxtf to %d", nxtf);
		}
                sprintf (syscmd, "mv DXDATO %s", rdbuf);
		free (rdbuf);
            } else {
		(void) system ("rm DXDADD*");
                sprintf (syscmd, "mv DXDATO DXDATM");
            }
            (void) system (syscmd);
	    /* get max and min in data */
	    /* average as for final plot */
	    i = 1;
	    a = 0.;
	    frequency = 0.;
	    m = 0;
	    df = (log10((double) fhigh) - log10((double) flow)) / smooth_bins;
	    logf = log10((double) nu[1]);
	    while (i < (max2 - 1)) {
		i++;
		while (log10((double) nu[i]) - logf > df) {
		    if (m != 0) {
			a /= m;
			frequency /= m;
			if (a * (1. - show_errors / (sigma * sqrt((double) m))) < botl) {
			    botl = a * (1. - show_errors / (sigma * sqrt((double) m)));
			}
			if (a * (1. + show_errors / (sigma * sqrt((double) m))) > upl)
			    upl = a * (1. + show_errors / (sigma * sqrt((double) m)));
			a = 0.;
			frequency = 0.;
			m = 0;
		    }
		    logf += df;
		}
		frequency += nu[i];
		if (iset == 0)
		    a += ps[i] * nu[i];
		else
		    a += psc[i] * nu[i];
		m++;
		if (smooth_bins == MAXSMOOTH) {
		    if (a * (1. - show_errors / sigma) < botl) {
			botl = a * (1. - show_errors / sigma);
		    }
		    if (a * (1. + show_errors / sigma) > upl)
			upl = a * (1. + show_errors / sigma);
		    frequency = a = 0.;
		    m = 0;
		}
	    }
	    if (m != 0) {
		a /= m;
		frequency /= m;
		if (a * (1. - show_errors / (sigma * sqrt((double) m))) < botl) {
		    botl = a * (1. - show_errors / (sigma * sqrt((double) m)));
		}
		if (a * (1. + show_errors / sigma) > upl)
		    upl = a * (1. + show_errors / sigma);
	    }
	    if (botl < 1.e-12 && 1.e-12 < upl)
		botl = 1.e-12;
	    a = log10((double) botl);
	    if (a < 0.)
		i = a - 1.;
	    else
		i = a;
	    botl = pow((double) 10., (double) i);
	    if (upl > 1.e6 && 1.e6 > botl)
		upl = 1.e6;
	    a = log10((double) upl);
	    if (a > 0)
		i = a + 1.;
	    else
		i = a;
	    upl = pow((double) 10., (double) i);
	    if (fplow_spec > 0.)
		botl = fplow_spec;
	    if (fphigh_spec > 0.)
		upl = fphigh_spec;

//	now generate pylab plot:
        FILE 	    *templatefile, *outfile, *oldoutfile;
// set up for python-generated plot
// using 3 files: template, OttoFFT.py, and data. 
        char    Template[120];
        char    PythonPlot[120];
        char    OldPythonPlot[120];
        char           *ch;
        char            pyvar[NSZ];

 	char    pyColor[10][12] = {"black", "blue", "darkgreen", "red", "cyan", "magenta", "darkorange", "brown", "violet", "lightblue"}; 

        strncpy (pyvar, variable, NSZ);
        if ((ch = memchr (pyvar, ' ', NSZ)) != NULL) {
	    printf (" return is %d\n", ch);
	    *ch = '\0';
        }
        strcpy(Template, (char *) getenv("XANADU"));
        strcat(Template, "/src/otto/OttoMEMTemplate.py");
        printf (" Template file is %s.\n", Template);
        templatefile = fopen(Template, "r");
//      sprintf (PythonPlot, "%sPlot.py", pyvar+1);
//      outfile = fopen(PythonPlot, "w");
        char *rdbuf;
        rdbuf = (char *) malloc(121 * sizeof(char));
	printf (" in otto, addmem is %d\n", addmem);
        if (addmem) {  // add to existing MEMPlot.py
            outfile = fopen ("MEMPlot2.py", "w");
            oldoutfile = fopen ("MEMPlot.py", "r");
            while (fgets (rdbuf, 120, oldoutfile) != NULL) {
		if (strncmp("## XX-----AAAAA-----XX", rdbuf, 22)) {
		    fputs (rdbuf, outfile);
                } else { break; }
            }		// add info for second plot line
            fprintf (outfile, "show_ps = 0\n");
            fprintf (outfile, "Var = np.empty ([3, DLEN], 'f')\n");
            fprintf (outfile, "VarName = [None]*3\n");
// 		skip first chunk in template file
            while (fgets(rdbuf, 120, templatefile) != NULL) {
                if (strncmp("XX-----BBBBB-----XX", rdbuf, 19)) {
                } else {
                    break;
                }
            }
//		then output second chunk
            while (fgets(rdbuf, 120, templatefile) != NULL) {
                if (strncmp("XX-----BBBBB-----XX", rdbuf, 19)) {
                    fputs (rdbuf, outfile);
                } else {
                    break;
                }
            }
	    sprintf (rdbuf, "DataFile = open ('./DXDADD%d', 'rb')\n", nxtf);
            fprintf (outfile, rdbuf);
            fprintf (outfile, "VarName[1] = \'nu\'\n");
            fprintf (outfile, "nu = np.empty (DLEN, 'f')\n");
            fprintf (outfile, "VarName[2] = \'ps\'\n");
            fprintf (outfile, "ps = np.empty (DLEN, 'f')\n");
            while (fgets(rdbuf, 120, templatefile) != NULL) {
                if (strncmp("XX-----BBBBB-----XX", rdbuf, 19)) {
                    fputs (rdbuf, outfile);
                } else {
                    break;
                }
            }
            fprintf (outfile, "    Vbuf = struct.unpack (3*'f', DataFile.read ((numvars)*4))\n");
            fprintf (outfile, "    Time[i] = Vbuf[0]\n");
            fprintf (outfile, "    nu[i] = Vbuf[1]\n");
            fprintf (outfile, "    ps[i] = Vbuf[2]\n");
            fprintf (outfile, "nu = np.ma.masked_where (nu == -32767., nu)\n");
            fprintf (outfile, "ps = np.ma.masked_where (ps == -32767., ps)\n");
            while (fgets(rdbuf, 120, templatefile) != NULL) {
                if (strncmp("XX-----BBBBB-----XX", rdbuf, 19)) {
                    fputs (rdbuf, outfile);
                } else {
                    break;
                }
            }
            while (fgets(rdbuf, 120, templatefile) != NULL) {
                if (strncmp("## XX-----CCCCC-----XX", rdbuf, 22)) {
                    fputs (rdbuf, outfile);
                } else {
                    break;
                }
            }
	    fprintf (outfile, "if show_smoothed:\n");
	    printf (" color index is %d\n", jcolor);
	    sprintf (rdbuf, "    pylab.plot (avenu, aveps, label='%s', color='%s', lw=2)\n", pyvar, pyColor[jcolor-1]);
    	    fprintf (outfile, rdbuf);
	    fprintf (outfile, "## XX-----AAAAA-----XX\n");
//		then add the rest from the original plot file:
            while (fgets (rdbuf, 120, oldoutfile) != NULL) {
		fputs (rdbuf, outfile);
            }
            free (rdbuf);
            printf (" ready to close outfile\n");
            (void) fclose(outfile);
	    sprintf (syscmd, "mv MEMPlot2.py MEMPlot.py");
	    (void) system (syscmd);
            // sprintf (syscmd, "ipython --pylab wx MEMPlot.py &"); 
            sprintf (syscmd, "python MEMPlot.py &"); 
            (void) system (syscmd);
        } else {
       	    // non-add run, remove addition data file:
	    if (temp_file = fopen ("DXDADD1", "r")) {
	        fclose (temp_file);
    	        sprintf (syscmd, "rm DXDADD*");
    	        (void) system (syscmd);
 	    }
            outfile = fopen ("MEMPlot.py", "w");
            while (fgets(rdbuf, 120, templatefile) != NULL) {
                if (strncmp("XX-----BBBBB-----XX", rdbuf, 19)) {
	            fputs (rdbuf, outfile);
                } else { break; }
            }
            fprintf (outfile, "numvars = 3\n");
            fprintf (outfile, "Date = %d\n", date);
            fprintf (outfile, "DLEN = %d\n", max2-1);
            fprintf (outfile, "Var = np.empty ([3, DLEN], 'f')\n");
            fprintf (outfile, "VarName = [None]*3\n");
            while (fgets(rdbuf, 120, templatefile) != NULL) {
                if (strncmp("XX-----BBBBB-----XX", rdbuf, 19)) {
                    fputs (rdbuf, outfile);
                } else {
                    break;
                }
            }
            fprintf (outfile, "DataFile = open('./DXDATM', 'rb')\n");
            fprintf (outfile, "VarName[1] = \'nu\'\n");
            fprintf (outfile, "nu = np.empty (DLEN, 'f')\n");
            fprintf (outfile, "VarName[2] = \'ps\'\n");
            fprintf (outfile, "ps = np.empty (DLEN, 'f')\n");
            while (fgets(rdbuf, 120, templatefile) != NULL) {
                if (strncmp("XX-----BBBBB-----XX", rdbuf, 19)) {
                    fputs (rdbuf, outfile);
                } else {
                    break;
                }
            }
            fprintf (outfile, "    Vbuf = struct.unpack (3*'f', DataFile.read ((numvars)*4))\n");
            fprintf (outfile, "    Time[i] = Vbuf[0]\n");
            fprintf (outfile, "    nu[i] = Vbuf[1]\n");
            fprintf (outfile, "    ps[i] = Vbuf[2]\n");
            fprintf (outfile, "nu = np.ma.masked_where (nu == -32767., nu)\n");
            fprintf (outfile, "ps = np.ma.masked_where (ps == -32767., ps)\n");
            while (fgets(rdbuf, 120, templatefile) != NULL) {
                if (strncmp("XX-----BBBBB-----XX", rdbuf, 19)) {
                    fputs (rdbuf, outfile);
                } else {
                    break;
                }
            }
            fprintf (outfile, "Fig = pylab.figure(figsize=(6.5,6))\n");
//      fprintf (outfile, "cid = Fig.canvas.mpl_connect ('button_press_event', onclick)\n");
            fprintf (outfile, "Fig.patch.set_facecolor ('#dddddd')\n");
            fprintf (outfile, "PanelL = Fig.add_subplot (111, axisbg='#ccddee')\n");
// 		// section for display of cursor values
            fprintf (outfile, "def format_coord2 (al, y):\n");
            fprintf (outfile, "    global DLEN, tas_average, avenu, aveps\n");
            fprintf (outfile, "    x = tas_average / al\n");
            fprintf (outfile, "    if x >= avenu[0] and x <= avenu[-1]:\n");
            fprintf (outfile, "        y1 = x\n");
            fprintf (outfile, "        y2 = al\n");
            fprintf (outfile, "        ix = 0\n");
            fprintf (outfile, "        while avenu[ix] < x and ix < len (avenu): ix += 1\n");
            fprintf (outfile, "        y3 = aveps[ix]\n");
            fprintf (outfile, "    else:\n");
            fprintf (outfile, "        y1 = -32767.\n");
            fprintf (outfile, "        y2 = -32767.\n");
            fprintf (outfile, "        y3 = -32767.\n");
            fprintf (outfile, "    return '(%%.2f, %%.2f, %%.2e)'%%(y1,y2,y3)\n");
            fprintf (outfile, "if show_ps:\n");
            fprintf (outfile, "    pylab.plot (nu, ps, label = 'PSD', color='red', lw=0.8)\n");
            fprintf (outfile, "    pylab.ylabel (r'%s:  f$\\ $P(f) = $-\\lambda P\\ (\\lambda)$', color='red', fontsize=20)\n", pyvar);
            fprintf (outfile, "tas_average=%f\n", tas_average);
            fprintf (outfile, "lateral = %d\n", spect_variable);
	    if (flow_spec > 0.) {
                fprintf (outfile, "flow=%f\n", flow_spec);
	    } else {
                fprintf (outfile, "flow=%f\n", flow);
	    }
            if (fhigh_spec > 0.) {
                fprintf (outfile, "fhigh=%f\n", fhigh_spec);
	    } else {
                fprintf (outfile, "fhigh=%f\n", fhigh);
	    }
	    fprintf (outfile, "fplow=%f\n", fplow_spec);
	    fprintf (outfile, "fphigh=%f\n", fphigh_spec);
            fprintf (outfile, "smooth_bins=%d\n", smooth_bins);
            fprintf (outfile, "variance = %.1f\n", variance);
            fprintf (outfile, "EDR = %.5f\n", eps_bar);
            while (fgets(rdbuf, 120, templatefile) != NULL) {
		//printf ("line is %s", rdbuf);
                if (strncmp("## XX-----CCCCC-----XX", rdbuf, 22)) {
                    fputs (rdbuf, outfile);
                } else {
                    break;
                }
            }
            fprintf (outfile, "if show_smoothed:\n");
            printf (" color index is %d\n", jcolor);
            sprintf (rdbuf, "    pylab.plot (avenu, aveps, label='%s', color='%s', lw=2)\n", pyvar, pyColor[jcolor-1]);
            fprintf (outfile, rdbuf);
            fgets(rdbuf, 120, templatefile);
            fgets(rdbuf, 120, templatefile);
            while (fgets(rdbuf, 120, templatefile) != NULL) {
                if (strncmp("XX-----BBBBB-----XX", rdbuf, 19)) {
                    fputs (rdbuf, outfile);
                } else {
                    break;
                }
            }

            fprintf (outfile, "PanelL.annotate ('Method: MEM; Variance=%.2f m$^2$s$^{-2}$, poles=%d\\nresolution=%.1e, bins for smoothing=%d\\nestimated eddy dissipation rate: %.1e m$^2$s$^{-3}$', xycoords='axes fraction', xy=(0.02, 0.05), backgroundcolor='lightyellow')\n", variance, poles, resolution, smooth_bins, eps_bar);
            fprintf (outfile, "PanelL.annotate ('%d %d--%d', xycoords='axes fraction', xy=(0.5,0.937), backgroundcolor='lightyellow')\n", date, start_time, end_time);
	    // need to fix this XXX
	    //sprintf (rdbuf, "%s: %s smoothed\\nred: %s", pyColor[icolor], pyvar);
	    sprintf (rdbuf, "PanelL.annotate ('%s', xycoords='axes fraction', xy=(0.6, 0.25), backgroundcolor='lightyellow')\n");
	    fprintf (outfile, rdbuf);
//	    fprintf (outfile, "PanelL.annotate ('blue: WIY smoothed\\nred: WIY', xycoords='axes fraction', xy=(0.6, 0.25), backgroundcolor='lightyellow')\n");
//          fprintf (outfile, "formatter = pylab.FuncFormatter(log_10_product)\n");
//          fprintf (outfile, "PanelL.yaxis.set_major_formatter(formatter)\n");
//          fprintf (outfile, "PanelL.xaxis.set_major_formatter(formatter)\n");
//          fprintf (outfile, "PanelL.yaxis.set_minor_formatter(formatter)\n");
//          fprintf (outfile, "if ymin > 0. and ymax > 0.: pylab.ylim ([ymin, ymax])\n");
//          fprintf (outfile, "ylim = PanelL.get_ylim()\n");
//          fprintf (outfile, "pylab.ylim(ylim)\n");
            free (rdbuf);
            fprintf (outfile, "pylab.grid()\n");
//            fprintf (outfile, "datacursor (draggable=True, formatter='{label}'.format)\n");
            fprintf (outfile, "flim = PanelL.get_xlim()\n");
            fprintf (outfile, "PanelT = PanelL.twiny()\n");
            fprintf (outfile, "PanelT.set_xlabel('wavelength [m]', color='black')\n");
            fprintf (outfile, "pylab.xlim ([tas_average / flim[0], tas_average / flim[1]])\n");
            fprintf (outfile, "PanelT.set_xscale ('log')\n");
            fprintf (outfile, "for tl in PanelT.get_xticklabels():\n");
            fprintf (outfile, "    tl.set_color('black')\n");
            fprintf (outfile, "PanelT.format_coord = format_coord2\n");
    
            fprintf (outfile, "pylab.savefig('MEMPlot.png', facecolor='#dddddd')\n");
            fprintf (outfile, "# pylab.show ()\n");
            printf (" ready to close outfile\n");
            (void) fclose(outfile);
            // sprintf (syscmd, "ipython --pylab wx MEMPlot.py &"); 
            sprintf (syscmd, "python MEMPlot.py &"); 
            (void) system (syscmd);
        }

	    ncar_(&flow, &fhigh, &botl, &upl, &ione, &ione, &ione, &ione, &ifour, &ione);
	    lablx_(":PGL:T:RU: :L:(:U:H:L:Z):U:", 27);
	    ilbl = 0;
	    if (iset == 0) {
	        if (!strncmp(variable, "WI", 2)) ilbl = 1;
	        if (!strncmp(variable, "LONW", 4)) ilbl = 1;
	        if (!strncmp(variable, "LATW", 4)) ilbl = 1;
	        if (!strncmp(variable, "LATW", 4)) ilbl = 1;
	    } else {
	        if (!strncmp(covariable, "WI", 2)) ilbl = 1;
	        if (!strncmp(covariable, "LONW", 4)) ilbl = 1;
	        if (!strncmp(covariable, "LATW", 4)) ilbl = 1;
	        if (!strncmp(covariable, "LATW", 4)) ilbl = 1;
	    }
	    if (ilbl)
	        lably_(":GL:T:RU:P(:GL:T:RU:) [m:S1:2 s:S2:-2  ]", 40);
	    else
	        lably_(":GL:T:RU:P(:GL:T:RU:)", 21);
	    add_lscale(&flow, &fhigh, &botl, &upl);
	    gcolor_(&icolor);
	    /* now average together values within smooth_bins */
	    for (isgn = -1; isgn <= 1; isgn++) {
		if (isgn && !show_errors)
		    continue;
		if (isgn) {
		    gcolor_(&ithree);
		    setusv_("LW", &thousand, 2);
		} else {
		    gcolor_(&icolor);
		    setusv_("LW", &two_thousand, 2);
		}
		df = (log10((double) fhigh) - log10((double) flow)) / smooth_bins;
		logf = log10((double) nu[1]);
		ifrst = 1;
		i = 0;
		a = 0.;
		frequency = 0.;
		m = 0;
		while (i < max2-1) {	// change 110805
		    i++;
		    while (log10((double) nu[i]) - logf > df) {
			if (m != 0) {
			    a /= m;
			    a *= (1. + isgn / (sigma * sqrt((double) m)));
			    frequency /= m;
			    if (frequency < flow) {ifrst = 1.;}
			    if (a < botl)
				a = botl;
			    if (ifrst) {
				ifrst = 0;
				frstpt_(&frequency, &a);
			    } else
				vector_(&frequency, &a);
			    a = 0.;
			    frequency = 0.;
			    m = 0;
			}
			logf += df;
		    }
		    frequency += nu[i];
		    if (iset == 0)
			a += ps[i] * nu[i];
		    else
			a += psc[i] * nu[i];
		    m++;
		    if (smooth_bins == MAXSMOOTH) {
			a *= (1. + isgn / sigma);
			if (ifrst) {
			    ifrst = 0;
			    frstpt_(&frequency, &a);
			} else
			    vector_(&frequency, &a);
			frequency = a = 0.;
			m = 0;
		    }
		}
		if (m != 0) {
		    a /= m;
		    a *= (1. + isgn / (sigma * sqrt((double) m)));
		    frequency /= m;
		    if (a < botl)
			a = botl;
		    vector_(&frequency, &a);
		}
	    }			/* end of isgn loop */
	    /* add constant-epsilon lines  */
	    gcolor_(&iseven);
	    setusv_("LW", &thousand, 2);
	    dashdb_(&dash_pat2);
	    if (iset == 0) {
		if (spect_variable == 0 || !strncmp (variable, "UX", 2))
		    ae = 0.15;
		else
		    ae = 0.2;
	    } else {
		if (cospect_variable == 0 || !strncmp (covariable, "UX", 2))
		    ae = 0.15;
		else
		    ae = 0.2;
	    }
	    for (i = -8; i < 0; i++) {
		if (i == -4) {
		    setusv_("LW", &two_thousand, 2);
		    dashdb_(&dash_pattern);
		} else if (i == -3) {
		    setusv_("LW", &thousand, 2);
		    dashdb_(&dash_pat2);
		}
		a = ae * pow((double) 10., (double) i * 0.66667);
		a *= pow(tas_average, 0.66667);
		a1 = a * pow((double) flow, -0.66667);
		if (a1 > upl) {
		    a1 = pow((double) 10., (double) i);
		    a1 *= tas_average;
		    a1 *= pow((double) upl / ae, (double) -1.5);
		    if (a1 < fhigh)
			frstd_(&a1, &upl);
		} else
		    frstd_(&flow, &a1);
		a2 = a * pow((double) fhigh, -0.66667);
		if (a2 < botl) {
		    a2 = pow((double) 10., (double) i);
		    a2 *= tas_average;
		    a2 *= pow((double) botl / ae, (double) -1.5);
		    if (a2 > flow && a1 < fhigh) {
			vectd_(&a2, &botl);
			if (i == -4) {
			    a1 = 1.20 * botl;
			    a = (log10((double) upl) - log10((double) botl))
				/ (log10((double) fhigh) - log10((double) flow));
			    a = -0.66667 / a;
			    a = atan((double) a) * 180. / 3.141592;
			    size = 8.;
			    cntr = 1.;
			    plchhq_(&a2, &a1,
			    "1:KGU:V:PRU:10:S2:-4                         ",
				    &size, &a, &cntr, 35);
			}
		    }
		} else {
		    if (a1 < fhigh) {
			vectd_(&fhigh, &a2);
			if (i == -4) {
			    a2 *= 1.20;
			    a = (log10((double) upl) - log10((double) botl))
				/ (log10((double) fhigh) - log10((double) flow));
			    a = -0.66667 / a;
			    a = atan((double) a) * 180. / 3.141592;
			    size = 8.;
			    cntr = 1.;
			    plchhq_(&fhigh, &a2,
			    "1:KGU:V:PRU:10:S2:-4                         ",
				    &size, &a, &cntr, 35);
			}
		    }
		}
	    }
	    setusv_("LW", &thousand, 2);
	    dashdb_(&dash_pat2);
	    for (i = -8; i <= 2; i++) {
		if (i == -4) {
		    setusv_("LW", &two_thousand, 2);
		    dashdb_(&dash_pattern);
		} else if (i == -3) {
		    setusv_("LW", &thousand, 2);
		    dashdb_(&dash_pat2);
		}
		a = pow((double) 10., (double) i);
		a1 = a * flow;
		if (a1 < botl) {
		    a1 = botl / a;
		    if (a1 > flow) {
			frstd_(&a1, &botl);
		    }
		} else
		    frstd_(&flow, &a1);
		a2 = a * fhigh;
		if (a2 > upl) {
		    a2 = upl / a;
		    if (a2 > flow)
			vectd_(&a2, &upl);
		} else {
		    if (a2 > botl)
			vectd_(&fhigh, &a2);
		}
	    }
	    setusv_("LW", &two_thousand, 2);
	    gcolor_(&ione);
	    a = fhigh * 0.9;
	    upl *= 0.5;
	    cntr = 1.;
	    size = 12.;
	    angd = 0.;
	    (void) sprintf(label, "VARIANCE=%12.4g", variance);
	    printf(" MEM total variance=%12.4g, <20-km=%12.4g\n", variance,
		   variance_20km);
	    label_length = 24;
	    /* plchhq_(&a, &upl, label, &size, &angd, &cntr, label_length); */
	    head_(&date, &start_time, &end_time, &dummy, 1);
	    titleline = 1;
	    if (iset == 0)
		fnote_(&titleline, variable, strlen(variable));
	    else
		fnote_(&titleline, covariable, strlen(covariable));
	    idline_(&titleline, label, strlen(label));
	    titleline = 2;
	    fnote_(&titleline, "MEM", 3);
	    (void) sprintf(label, "poles=%5d, resln=%8.6f, smooth bins=%5d",
			   poles, resolution, smooth_bins);
	    idline_(&titleline, label, 46);
	    titleline = 3;
	    if (iset == 0) {
		if (!strncmp(variable, "WIND", 4)) {
		    sprintf(label, "wind in %d degree direction", proj_degrees);
		    idline_(&titleline, label, strlen(label));
		} else {
		    sprintf(label, ":GL:E:RU:=%.3e :GL:,:RU:%.3e", eps_bar, eps2_bar);
		    idline_(&titleline, label, strlen(label));
		}
	    } else {
		if (!strncmp(covariable, "WIND", 4)) {
		    sprintf(label, "wind in %d degree direction", cproj_degrees);
		    idline_(&titleline, label, strlen(label));
		} else {
		    sprintf(label, ":GL:E:RU:=%.3e :GL:,:RU:%.3e", eps_bar, eps2_bar);
		    idline_(&titleline, label, strlen(label));
		}
	    }
	    gcolor_(&ione);

	    /*
	     * this strange call forces color reset, needed if plot frames
	     * are overlaid by 'med'
	     */
	    line_(&thousand, &thousand, &thousand, &thousand);
	    frame_();
	}
	printf(" after fP(f) section\n");
	if (((ishowmem & 8) && !iset) || ((ishowcmem & 8) && iset)) {
	    /* set limits */
	    upl = 1.e-10;
	    botl = 1.e10;
	    /* set limits from extremes in data: */
	    if (iset == 0) {
		if (spect_variable == 0  || !strncmp (variable, "UX", 2))
		    ae = 2.0;
		else
		    ae = 1.5;
	    } else {
		if (cospect_variable == 0 || !strncmp (covariable, "UX", 2))
		    ae = 2.0;
		else
		    ae = 1.5;
	    }
	    ae = pow(ae, (double) 1.5) * TWOPI / tas_average;
	    eps_bar = 0.;
	    eps2_bar = 0.;
	    m_eps = 0;
	    for (i = 1; i <  max2; i++) {	/* convert to equiv of eddy
						 * diss. rate */
//		if (nu[i] >= nu[max2-1] / 10.) {
  		if ((nu[i] > 1.0) && (nu[i] <= 8.)) {
		    if (iset == 0) {
//			a = pow((double) ps[i], (double) (1.5)) * ae
//			    * pow((double) nu[i], 2.5);
			a = ps[i] * pow ((double) nu[i], FIVETHIRDS);  // corr Oct 2013
		    } else {
//			a = pow((double) psc[i], (double) (1.5)) * ae
//			    * pow((double) nu[i], 2.5);
			a = psc[i] * pow ((double) nu[i], FIVETHIRDS);  // corr Oct 2013
		    }
		    eps_bar += a;
		    eps2_bar += a * a;
		    m_eps++;
		}
		if (iset == 0) {
//		    b = pow((double) ps[i], (double) (1.5)) * ae
//			* pow((double) nu[i], 2.5);
		    b = ps[i] * pow ((double) nu[i], FIVETHIRDS);
		} else {
//		    b = pow((double) psc[i], (double) (1.5)) * ae
//			* pow((double) nu[i], 2.5);
		    b = psc[i] * pow ((double) nu[i], FIVETHIRDS);
		}
		b = ae * pow ((double) b, (double) 1.5);  // Oct 2013
		if (botl > b * (1. - show_errors / sigma))
		    botl = b * (1. - show_errors / sigma);
		if (upl < b * (1. + show_errors / sigma))
		    upl = b * (1. + show_errors / sigma);
	    }
	    if (m_eps > 0) {
		eps_bar /= m_eps;
		eps2_bar /= m_eps;
		eps2_bar = sqrt((eps2_bar - eps_bar * eps_bar));
		eps_bar = ae * pow ((double) eps_bar, 1.5); // corr Oct 2013
		eps2_bar = ae * pow ((double) eps2_bar, 1.5);  // corr Oct 2013
		printf(" MEM average eddy dissipation rate=%.5e +/- %.5e\n",
		       eps_bar, eps2_bar);
	    }
	    if (botl < 1.e-12 && 1.e-12 < upl)
		botl = 1.e-12;
	    a = log10((double) botl);
	    if (a < 0)
		i = a - 1.;
	    else
		i = a;
	    botl = pow((double) 10., (double) i);
	    if (upl > 1.e6 && 1.e6 > botl)
		upl = 1.e6;
	    a = log10((double) upl);
	    if (a > 0)
		i = a + 1.;
	    else
		i = a;
	    upl = pow((double) 10., (double) i);
	    if (epslow_spec > 0.)
		botl = epslow_spec;
	    if (epshigh_spec > 0.)
		upl = epshigh_spec;
	    ncar_(&flow, &fhigh, &botl, &upl, &ione, &ione, &ione, &ione, &ifour, &ione);
	    lablx_(":PGL:T:RU: :L:(:U:H:L:Z):U:", 27);
	    lably_("E(:GL:T:RU:) [m:S1:2 s:S2:-3  ]", 31);
	    add_lscale(&flow, &fhigh, &botl, &upl);
	    gcolor_(&icolor);
	    /* now average together values within smooth_bins */
	    for (isgn = -1; isgn <= 1; isgn++) {
		if (isgn && !show_errors)
		    continue;
		if (isgn) {
		    gcolor_(&ithree);
		    setusv_("LW", &thousand, 2);
		} else {
		    gcolor_(&icolor);
		    setusv_("LW", &two_thousand, 2);
		}
		df = (log10((double) fhigh) - log10((double) flow)) / smooth_bins;
		logf = log10((double) nu[1]);
		ifrst = 1;
		i = 0;
		a = 0.;
		frequency = 0.;
		m = 0;
		while (i < max2) {
		    i++;
		    while (log10((double) nu[i]) - logf > df) {
			if (m != 0) {
			    a /= m;
			    a *= (1. + 1.5 * isgn / (sigma * sqrt((double) m)));
			    a = ae * pow (a, 1.5);
			    frequency /= m;
			    if (frequency < flow) {ifrst = 1.;}
			    if (a > upl)
				a = upl;
			    if (a < botl)
				a = botl;
			    if (ifrst) {
				ifrst = 0;
				frstpt_(&frequency, &a);
			    } else
				vector_(&frequency, &a);
			    a = 0.;
			    frequency = 0.;
			    m = 0;
			}
			logf += df;
		    }
		    frequency += nu[i];
		    if (iset == 0)
//			a = pow((double) ps[i], (double) (1.5)) * ae
//			    * pow((double) nu[i], 2.5);
				// no correction needed, indiv pts, Oct 2013
				// but made for consistency:
			a += ps[i] * pow ((double) nu[i], (double) FIVETHIRDS);
		    else
//			a = pow((double) psc[i], (double) (1.5)) * ae
//			    * pow((double) nu[i], 2.5);
			a = psc[i] * pow ((double) nu[i], (double) FIVETHIRDS);
		    m++;
		    if (smooth_bins == MAXSMOOTH) {
			a *= (1. + 1.5 * isgn / sigma);
			a = ae * pow (a, 1.5);
			if (a > upl)
			    a = upl;
			if (a < botl)
			    a = botl;
			if (ifrst) {
			    ifrst = 0;
			    frstpt_(&frequency, &a);
			} else
			    vector_(&frequency, &a);
			frequency = a = 0.;
			m = 0;
		    }
		}
		if (m != 0) {
		    a /= m;
		    frequency /= m;
		    a *= (1. + 1.5 * isgn / (sigma * sqrt((double) m)));
		    a = ae * pow (a, 1.5);
		    if (a > upl)
			a = upl;
		    if (a < botl)
			a = botl;
		    vector_(&frequency, &a);
		}
	    }			/* end of isgn loop */
	    setusv_("LW", &two_thousand, 2);
	    gcolor_(&iseven);
	    frstpt_(nu + max2-1, &eps_bar);
	    a = nu[max2-1] / 10.;
	    vector_(&a, &eps_bar);
	    gcolor_(&ione);
	    a = fhigh * 0.9;
	    upl *= 0.5;
	    cntr = 1.;
	    size = 12.;
	    angd = 0.;
//	    (void) sprintf(label, "VARIANCE=%15.7g", variance);
  	    (void) sprintf(label, "EDR=%15.7g", eps_bar);
	    label_length = 24;
	    /* plchhq_(&a, &upl, label, &size, &angd, &cntr, label_length); */
	    head_(&date, &start_time, &end_time, &dummy, 1);
	    titleline = 1;
	    if (iset == 0)
		fnote_(&titleline, variable, strlen(variable));
	    else
		fnote_(&titleline, covariable, strlen(covariable));
	    idline_(&titleline, label, strlen(label));
	    titleline = 2;
	    fnote_(&titleline, "MEM", 3);
	    sprintf(label, ":GL:E:RU:=%.3e :GL:,:RU:%.3e", eps_bar, eps2_bar);
	    idline_(&titleline, label, strlen(label));
	    if (iset == 0) {
		if (!strncmp(variable, "WIND", 4)) {
		    titleline = 3;
		    sprintf(label, "wind in %d degree direction", proj_degrees);
		    idline_(&titleline, label, strlen(label));
		}
	    } else {
		if (!strncmp(covariable, "WIND", 4)) {
		    titleline = 3;
		    sprintf(label, "wind in %d degree direction", cproj_degrees);
		    idline_(&titleline, label, strlen(label));
		}
	    }
	    gcolor_(&ione);

	    /*
	     * this strange call forces color reset, needed if plot frames
	     * are overlaid by 'med'
	     */
	    line_(&thousand, &thousand, &thousand, &thousand);
	    frame_();
	}
	printf(" after epsilon-dependent plot\n");
	if ((ishowmem & 16) && !iset) {	/* forecast plot, only for "variable" */
	    fbuffer_size = buffer_size * 1.2;
	    if (fcst_buffer == NULL) {
    printf(" attempt #5 to allocate %d float words\n", fbuffer_size);
                fcst_buffer = (float *) GetMemory((size_t) fbuffer_size * sizeof(float));
	    }
	    fmin = 1. / points;
	    fmax = 0.5;
	    fmin = log((double) fmin);
	    fmax = log((double) fmax);
	    df = (fmax - fmin) / max2;
    	    t = 0.;
    	    data_index = buffer_size ;
#if(STOCKS)
			/* reverse order for normal forecast appearance */
    	    while (data_index > 0) {
	        *(fcst_buffer + data_index - 1) = *(data_buffer + buffer_size - data_index);
		data_index--;
	    }
#else
	    while (data_index > 0) {
		*(fcst_buffer + data_index - 1) = *(data_buffer + data_index - 1);
		data_index--;
	    }
#endif
		/* now construct forecast 20% into the future: */
#if(ITERATING)
#else
            fprintf(fcstfile, "  FORECAST \n\n");
	    fprintf(fcstfile, " RMS for linear prediction coefficients: %f\n", pm);
#endif
  	    value = *(fcst_buffer + data_size -1);
#if(ITERATING)
#else
  	    value += (mean + trend * (duration/2.));
#endif
#if(STOCKS)
	    value -= trend * duration;
	    fprintf(fcstfile, " last value in time series is %f\n", value);
#endif
#if(ITERATING)
	    fprintf(fcstfile, "%f", value);
#endif
	    printf(" last value in time series is %f; data_size=%d\n", value,data_size);
	    value_Last = pow(10., value);
#if(STOCKS)
	    ptrend = -1. * trend;
#else
	    ptrend = trend;
#endif
#if(STOCKS)
	    fprintf(fcstfile, "mean = %f, trend = %f\n", mean, ptrend);
	    if(variable[0] == 'L') {
	        value = (pow(10., 250.*ptrend)-1.)*100.;
		localRMS = (pow(10., LogRMS) - 1.) * 100.;
	        fprintf(fcstfile, "annual growth (linear trend)=%3.0f\% RMS=%3.0f\%\n", value, localRMS);
	    }
#endif
#if(ITERATING)
	    for (j = 0; j < 10; j++) {
#else
  	    for (j = 0; j < (fbuffer_size-buffer_size); j++) {
#endif

		value = 0.;
	        for (i = 1; i <= poles; i++) {
	            value += cof[i] * *(fcst_buffer + data_size + j - i);
		}
		*(fcst_buffer + data_size + j) = value;
	        value += (mean + ptrend * (duration/2. + (j+1)));
#if(STOCKS)
  		fprintf(fcstfile, " forecast %d day = %f\n", j, value);
#endif
#if(ITERATING)
  		fprintf(fcstfile, ";%f", value);
#endif
		if (j == pday1) {
		    value_Day1 = pow(10., value);
//		    sprintf(labelF1, "Forecast %6.2f after %d days", value, j);
		    printf("%s Forecast %6.2f after %d days, ", (variable+1), value, j);
		}
		if (j == pday2) {
		    value_Day2 = pow(10., value);
		    sprintf(labelF2, "Forecast %6.2f after %d days, %6.2f after %d days", value_Day1, pday1, value_Day2, j);
		    printf(" %6.2f after %d days\n", value, j);
		}
#if(ITERATING)
		if (j < 10) {printf(" %6.2f after %d days\n", value, j+1);}
#endif
	    }
#if(ITERATING)
	    fprintf(fcstfile, "\n");
#else
	    plotfcst(fcst_buffer, mean, ptrend, variable, proj_degrees, &vmin, &vmax);
#endif
	}
#if(ITERATING)
    data_size--;
    } while (--iterator > 0);

    data_size += 1000/dhd1.ideltt;
#endif
    }
    printf(" close fcstfile, return is %d\n", fclose(fcstfile));
//  (void) fclose(fcstfile);
    printf(" ready for cospectrum and quadrature...\n");
    /* now plot cospectrum and quadrature */
    if (ishowcmem & 48) {
					/* get autoregressive estimates */

	int ijij;
	for (ijij = 0; ijij<1; ijij++) {  // this is just for "continue" below
	cofcv = vector(1, 4 * poles);
        mcar_(&itwo, &data_size, &poles, data_buffer, co_data_buffer,
	    pmcv, cofcv, &status_mcar);
	if (status_mcar == 4) {continue;}	// branch for too many points
	cospectrum = vector(1, max2);	
	quadrature = vector(1, max2);	
	fmin = 1. / points;
	fmax = 0.5;
	fmin = log((double) fmin);
	fmax = log((double) fmax);
	df = (fmax - fmin) / max2;
	fdt = 0.0;
	for (i = 1; i <= max2; i++) {
	    last_f = fdt;
	    fdt = fmin + df * i;
	    fdt = exp((double) fdt);
 	    mcarpsd_(&poles, &fdt, pmcv, cofcv, &ps[i], &psc[i], &cospectrum[i], &quadrature[i],
		&status_mcar);
	    cospectrum[i] *= 2. / sample_frequency;
	    quadrature[i] *= 2. / sample_frequency;
	    ps[i]  *= 2. / sample_frequency;
	    psc[i] *= 2. / sample_frequency;
#if(0)
	    printf(" nu=%f, ps=%f, psc=%f\n", nu[i], ps[i], psc[i]);
	    printf(" cos=%f, q=%f\n",
		   cospectrum[i], quadrature[i]);
	    a = sqrt((double) ((cospectrum[i] * cospectrum[i] + quadrature[i] * quadrature[i]) / (ps[i] * psc[i])));
	    printf(" coherence=%f, phase", a);
	    a = atan2(quadrature[i], cospectrum[i]) * 180. / PI;
	    printf(" %f\n", a);
#endif
	}
	for (iset = 0; iset < 2; iset++) {
	    if (!(ishowcmem & 16)) {continue;}
	    upl = -1.e10, botl = 1.e10;
	    flux = 0.;
	    fdt = 0.;
	    df = (fmax - fmin) / max2;
	    for (i = 1; i <= max2; i++) {
		last_f = fdt;
		fdt = fmin + df * i;
		fdt = exp((double) fdt);
		nu[i] = fdt * sample_frequency;
		if (nu[i] > flux_min_freq && nu[i] <= flux_max_freq)
		    flux += cospectrum[i] * (fdt - last_f) * sample_frequency;
	    }

	    /* now average together values within smooth_bins */
	    df = (log10((double) fhigh) - log10((double) flow)) / smooth_bins;
	    for (isgn = -1; isgn <= 1; isgn++) {
		if (isgn && !show_errors)
		    continue;
		logf = log10((double) nu[1]);
		i = 0;
		a = 0.;
		m = 0;
		while (i++ < max2) {
		    while (log10((double) nu[i]) - logf > df) {
			if (m) {
			    a /= m;
			    if (a * (1. + isgn / (sigma * sqrt((double) m))) < botl) {
				botl = a * (1. + isgn / (sigma * sqrt((double) m)));
			    }
			    if (a * (1. + isgn / (sigma * sqrt((double) m))) > upl) {
				upl = a * (1. + isgn / (sigma * sqrt((double) m)));
			    }
			    a = 0.;
			    m = 0;
			}
			logf += df;
		    }
		    if (iset == 0)
			a += nu[i] * cospectrum[i];
		    else
			a += nu[i] * quadrature[i];
		    m++;
		    if (smooth_bins == MAXSMOOTH) {
			if (a * (1. + isgn / sigma) < botl) {
			    botl = a * (1. + isgn / sigma);
			}
			if (a > upl * (1. + isgn / sigma)) {
			    upl = a * (1. + isgn / sigma);
			}
			a = 0.;
			m = 0;
		    }
		}
		if (m) {
		    a /= m;
		    if (a * (1. + isgn / (sigma * sqrt((double) m))) < botl) {
			botl = a * (1. + isgn / (sigma * sqrt((double) m)));
		    }
		    if (a * (1. + isgn / (sigma * sqrt((double) m))) > upl) {
			upl = a * (1. + isgn / (sigma * sqrt((double) m)));
		    }
		}
		ipc = i * 10 / max2 + 60;
		if (ipc != ipcold) {
		    ipcold = ipc;
		}
	    }
	    flow = 0.001;
	    fhigh = 10.;
	    if (sample_frequency > 10.) {
		fhigh *= 10.;
	    }
	    if (flow_spec > 0.)
		flow = flow_spec;
	    if (fhigh_spec > 0.)
		fhigh = fhigh_spec;

	    set_limits(&botl, &upl);
	    if (cplow_spec != 0.)
		botl = cplow_spec;
	    if (cphigh_spec != 0.)
		upl = cphigh_spec;
	    ncar_(&flow, &fhigh, &botl, &upl, &ione, &ione, &ifour, &itwo, &ithree, &ione);
	    lablx_(":PGL:T:RU: :L:(:U:H:L:Z):U:", 27);
/*	lably_("P(:GL:T:RU:) [m:S1:2 s:S2:-1  ]", 31); */
	    ilbl = 0;
	    if (!strncmp(covariable, "WI", 2)) ilbl++;
	    if (!strncmp(covariable, "LONW", 4)) ilbl++;
	    if (!strncmp(covariable, "LATW", 4)) ilbl++;
	    if (!strncmp(variable, "WI", 2)) ilbl++;
	    if (!strncmp(variable, "LONW", 4)) ilbl++;
	    if (!strncmp(variable, "LATW", 4)) ilbl++;
	    if (ilbl < 2) ilbl = 0;
	    if (ilbl) {
	        if (iset == 0)
		    lably_(":GL:T:RU:C(:GL:T:RU: [m:S1:2 s:S2:-2  :]", 40);
	        else
		    lably_(":GL:T:RU:Q(:GL:T:RU: [m:S1:2 s:S2:-2:  ]", 40);
	    } else {
	        if (iset == 0)
		    lably_(":GL:T:RU:C(:GL:T:RU:)          ", 31);
	        else
		    lably_(":GL:T:RU:Q(:GL:T:RU:)          ", 31);
	    }
	    add_lscale(&flow, &fhigh, &botl, &upl);
	    df = (log10((double) fhigh) - log10((double) flow)) / smooth_bins;
	    for (isgn = -1; isgn <= 1; isgn++) {
		if (isgn && !show_errors)
		    continue;
		if (isgn) {
		    gcolor_(&ithree);
		    setusv_("LW", &thousand, 2);
		} else {
		    gcolor_(&icolor);
		    setusv_("LW", &two_thousand, 2);
		}
		/* now average together values within smooth_bins */
		logf = log10((double) nu[1]);
		ifrst = 1;
		i = 1;
		a = 0.;
		frequency = 0.;
		m = 0;
		while (i++ < max2) {
		    while (log10((double) nu[i]) - logf > df) {
			if (m) {
			    frequency /= m;
			    a *= (1. + isgn / (sigma * sqrt((double) m))) / m;
			    if (a < botl)
				a = botl;
			    if (ifrst) {
				ifrst = 0;
				frstpt_(&frequency, &a);
			    } else
				vector_(&frequency, &a);
			    a = 0.;
			    frequency = 0.;
			    m = 0;
			}
			logf += df;
		    }
		    frequency += nu[i];
		    if (iset == 0)
			a += nu[i] * cospectrum[i];
		    else
			a += nu[i] * quadrature[i];
		    m++;
		    if (smooth_bins == MAXSMOOTH) {
			a *= (1. + isgn / sigma);
			if (ifrst) {
			    ifrst = 0;
			    frstpt_(&frequency, &a);
			} else
			    vector_(&frequency, &a);
			frequency = 0.;
			a = 0.;
			m = 0;
		    }
		}
		if (m) {
		    a *= (1. + isgn / (sigma * sqrt((double) m))) / m;
		    frequency /= m;
		    if (a < botl)
			a = botl;
		    vector_(&frequency, &a);
		}
	    }
	    setusv_("LW", &two_thousand, 2);
	    gcolor_(&ione);
	    a = fhigh * 0.9;
	    printf("FLUX=%15.7g, limits %f - %f Hz\n", 
		flux, flux_min_freq, flux_max_freq);
	    (void) sprintf(label, "FLUX=%15.7g, limits %f - %f Hz", 
		flux, flux_min_freq, flux_max_freq);
	    titleline = 1;
	    idline_(&titleline, label, strlen(label));
	    head_(&date, &start_time, &end_time, &dummy, 1);
	    sprintf(label, "%s x %s", variable, covariable);
	    fnote_(&titleline, label, strlen(label));
	    titleline = 2;
	    sprintf(label, "MEM, %s delayed %d ms", variable, delay);
	    fnote_(&titleline, label, strlen(label));
	    (void) sprintf(label, "poles=%5d, resln=%8.6f, smooth bins=%5d",
			   poles, resolution, smooth_bins);
	    idline_(&titleline, label, 46);
	    if (iset == 0) {
		if (!strncmp(variable, "WIND", 4)) {
		    titleline = 3;
		    sprintf(label, "wind in %d degree direction", proj_degrees);
		    idline_(&titleline, label, strlen(label));
		}
	    } else {
		if (!strncmp(covariable, "WIND", 4)) {
		    titleline = 3;
		    sprintf(label, "wind in %d degree direction", cproj_degrees);
		    idline_(&titleline, label, strlen(label));
		}
	    }
	    if (!strncmp(variable, "WIND", 4)) {
		titleline = 3;
		sprintf(label, "wind in %d degree direction", proj_degrees);
		idline_(&titleline, label, strlen(label));
	    }
	    if (!strncmp(covariable, "WIND", 4)) {
		titleline = 3;
		sprintf(label, "wind in %d degree direction", cproj_degrees);
		idline_(&titleline, label, strlen(label));
	    }
	    gcolor_(&ione);

	    /*
	     * this strange call forces color reset, needed if plot frames
	     * are overlaid by 'med'
	     */
	    line_(&thousand, &thousand, &thousand, &thousand);
	    frame_();
	}
	}
    }
    printf(" ready to call cphase, status_mcar=%d\n", status_mcar);
    if ((ishowcmem & 32) && !status_mcar) {
	variance_factor = sigma;
	cphase(nu, cospectrum, quadrature, ps, psc, variance_factor,
	       spect_variable, variable, cospect_variable, covariable,
	       proj_degrees, cproj_degrees);
    }
    free(data_buffer);
    free(fcst_buffer);
    free_vector(cof, 1, poles);
    free_vector(ps, 1, max2);
    free_vector(nu, 1, max2);
    if (ishowp) {free(cc_data_buffer);}
    if (ishowcmem) {
	free(co_data_buffer);
	free_vector(cofc, 1, poles);
	free_vector(psc, 1, max2);
	if ((ishowcmem & 16) && !(status_mcar)) {
	    free_vector(cospectrum, 1, max2);
	    free_vector(quadrature, 1, max2);
	}
    }
    return;
}

void
getdata(pmode)
    int             pmode;	/* 1 to pad data to power-of-2 size */

{


    int             i, j, k, nv;
    int             ifl = 0;
    int             jfl, ns;
    int             stime;
    int             mode, cmode;
    int             nbytes, ierr;
    int             xran = 33651;
    int             nloop, ipc, ipcold = 0;
    static int	    first_ran = TRUE;

    float	    time_stop;
    float           wx, wy, vx, vy;
    float	    midpoint;
    float	    offset, co_offset;
    double          heading = -999., wd, proj, cproj;
    int             missing_data;


    if (first_ran) {
        srand(xran);
	first_ran = FALSE;
    }
		// size needs +1 increment to include both end points
		// change 160211: by convention, exclude end pt
    data_size = (int) (duration / (0.001 * dhd1.ideltt));
    if (pmode == 1) {
	i = data_size / segment_length + 1;
	buffer_size = i * segment_length + segment_length / 2;
    } else
	buffer_size = data_size;
    printf(" attempt #6 to allocate %d float words\n", buffer_size);
    printf(" data_size=%d, duration=%f, segment_length=%d\n", data_size, duration, segment_length);
    data_buffer = (float *) GetMemory((size_t) buffer_size * sizeof(float));
    stime = -start_time;
    nv = 0;
    points = 0;
    mean = 0.;
    trend = 0.;
    filtered_5by3 = 0.;
    filtered_value = 0.;
    if (!strncmp(variable, "LONW", 4)) {
	mode = 1;
    } else if (!strncmp(variable, "LATW", 4)) {
	mode = 2;
    } else if (!strncmp(variable, "WIND", 4)) {
	mode = 3;
/*	printf(" enter direction for wind projection (deg):");
 *	scanf("%f",&projection);
 */
	proj = proj_degrees * CRADEG;
    } else if (!strncmp(variable, "SIM ", 4)) {
	mode = 4;
    } else {
	strncpy(arcbuf[4], variable, NSZ);
	mode = 0;
    }
    if (ishowcfft || ishowcmem || ishowp) {
	co_mean = 0.;
	co_trend = 0.;
    printf(" attempt #7 to allocate %d float words\n", buffer_size);
	co_data_buffer = (float *) GetMemory((size_t) buffer_size * sizeof(float));
	if (ishowp) {
    printf(" attempt #8 to allocate %d float words\n", buffer_size);
	    cc_data_buffer = (float *) GetMemory((size_t) buffer_size * sizeof(float));
	    strncpy(covariable, variable, NSZ);
	    covariable[0] = 'N';
	    strncpy(ccvariable, variable, NSZ);
	    ccvariable[0] = 'C';
	}
	if (!strncmp(covariable, "LONW", 4)) {
	    cmode = 1;
	} else if (!strncmp(covariable, "LATW", 4)) {
	    cmode = 2;
	} else if (!strncmp(covariable, "WIND", 4)) {
	    cmode = 3;
	    cproj = cproj_degrees * CRADEG;
	} else if (!strncmp(covariable, "SIM ", 4)) {
	    cmode = 4;
	} else if (!strncmp(covariable, "NONE", 4)) {
	    cmode = 5;
	} else {
	    strncpy(arcbuf[5], covariable, NSZ);
	    cmode = 0;
	}
	cvmin = 1.e30;
	cvmax = -1.e30;
    }
    vmin = 1.e30;
    vmax = -1.e30;
    tas_average = 0.;
    nloop = 0;
    time_stop = fsect( ftsec ((float) end_time) + 0.001 * (delay + 1)) + 0.01;
    dpoints = 0;
    while (1) {
	jfl = arc_(filename, &nv, buf, &stime, arcbuf);
  	if (points == 0) {
	    printf(" first value from ARC is %f\n", buf[6]);
	    printf(" time_stop=%f jfl=%d\n", time_stop, jfl);
	}
	if (date == 0)
	    date = (int) buf[0];
	if (jfl < 0)
	    break;
 	if (buf[1] > time_stop)
	    break;
	if (!(++nloop % 100)) {	/* note, only 0--20% displayed here */
	    ipc = (itsec((int) buf[1]) - itsec(start_time)) * 20 /
		(itsec(end_time) - itsec(start_time));
	    if (ipc != ipcold) {
		ipcold = ipc;
	    }
	}
	missing_data = FALSE;
	if (mode == 1) {
	    if (heading == -999.) {
		heading = buf[4] * CRADEG;
		if (buf[4] == MISSING_DATA)
		    missing_data = TRUE;
	    }
	    vx = sin(heading);
	    vy = cos(heading);
	    wd = buf[2] * CRADEG;
	    wx = sin(wd);
	    wy = cos(wd);
	    value = -buf[3] * (wx * vx + wy * vy);
	    if (buf[2] == MISSING_DATA || buf[3] == MISSING_DATA)
		missing_data = TRUE;
	} else if (mode == 2) {
	    if (buf[2] == MISSING_DATA || buf[3] == MISSING_DATA || buf[4] == MISSING_DATA)
		missing_data = TRUE;
	    heading = buf[4] * CRADEG;
	    vx = sin(heading);
	    vy = cos(heading);
	    wd = buf[2] * CRADEG;
	    wx = sin(wd);
	    wy = cos(wd);
	    value = buf[3] * (wx * vy - wy * vx);
	} else if (mode == 3) {
	    vx = sin(proj);
	    vy = cos(proj);
	    wd = buf[2] * CRADEG;
	    if (buf[2] == MISSING_DATA || buf[3] == MISSING_DATA)
		missing_data = TRUE;
	    wx = sin(wd);
	    wy = cos(wd);
	    value = -buf[3] * (wx * vx + wy * vy);
	} else if (mode == 4) {
	    filtered_5by3 += (((float) rand() / RAND_MAX-0.5) * 1000. * random_5by3 - filtered_5by3) / (sample_frequency * 100.);
	    value = filtered_5by3 + ((float) rand() / RAND_MAX-0.5) * random_noise;	
	    value += cosine_amplitude * cos((double) points * TWOPI / (cosine_period * sample_frequency));
		// special section to filter the results
	    if (filter_tau > 0.) {
	        filtered_value += (value - filtered_value) / (sample_frequency * filter_tau) ;
	        value = filtered_value;
	    }
	} else {
	    value = buf[6];
	    if (value == MISSING_DATA)
		missing_data = TRUE;
	}
	if (ishowcmem || ishowp || ishowcfft) {
	    if (cmode == 1) {
		if (heading == -999.) {
		    heading = buf[4] * CRADEG;
		    if (buf[4] == MISSING_DATA)
			missing_data = TRUE;
		}
		vx = sin(heading);
		vy = cos(heading);
		wd = buf[2] * CRADEG;
		wx = sin(wd);
		wy = cos(wd);
		cvalue = -buf[3] * (wx * vx + wy * vy);
	        if (buf[2] == MISSING_DATA || buf[3] == MISSING_DATA)
		    missing_data = TRUE;
	    } else if (cmode == 2) {
		heading = buf[4] * CRADEG;
		vx = sin(heading);
		vy = cos(heading);
		wd = buf[2] * CRADEG;
		wx = sin(wd);
		wy = cos(wd);
		cvalue = buf[3] * (wx * vy - wy * vx);
	        if (buf[2] == MISSING_DATA || buf[3] == MISSING_DATA || buf[4] == MISSING_DATA)
		    missing_data = TRUE;
	    } else if (cmode == 3) {
		vx = sin(cproj);
		vy = cos(cproj);
		wd = buf[2] * CRADEG;
		wx = sin(wd);
		wy = cos(wd);
		cvalue = -buf[3] * (wx * vx + wy * vy);
	        if (buf[2] == MISSING_DATA || buf[3] == MISSING_DATA)
		    missing_data = TRUE;
	    } else if (cmode == 4) {
		cvalue = (rand() / RAND_MAX - 0.5) * co_random_noise;
		cvalue += co_cosine_amplitude * cos((double) (points * TWOPI / (sample_frequency * co_cosine_period) + co_cosine_phase * TWOPI / 360.));
	    } else if (cmode == 5) {	// "None" entry -- skip
	    } else {
		cvalue = buf[7];
//		printf(" cvalue is %f\n", cvalue);
	        if (cvalue == MISSING_DATA)
		    missing_data = TRUE;
	    }
	}
	if (missing_data) {
	    value = MISSING_DATA;
	    cvalue = MISSING_DATA;
	} else
	    dpoints++; 
	/* note time=tstart+points*dhd1.ideltt/1000 */
	/* need protection or correction for time gaps */
//      midpoint = ((duration - 1.) / (0.002 * dhd1.ideltt));
        midpoint = ((duration) / (0.002 * dhd1.ideltt));
	if (points < data_size) {
	    if (dpoints == 1) {
		offset = value;
	    }
	    *(data_buffer + points) = value;	/* load into storage */
	    if (value != MISSING_DATA) {
		mean += (value - offset);
	        trend += (value - offset) * (points - midpoint);
	        if (value > vmax)
	            vmax = value;
	        if (value < vmin)
	            vmin = value;
	        tas_average += buf[5];
	    }
	}
	if (ishowcmem || ishowp || ishowcfft) {
	    if (dpoints == 1) {
		co_offset = cvalue;
	    }
	    if (points - idelay >= 0) {
	        *(co_data_buffer + points - idelay) = cvalue;
	        if (ishowp) {*(cc_data_buffer + points - idelay) = buf[8];}
		if (cvalue != MISSING_DATA) {
	            co_mean += (cvalue - co_offset);
	            co_trend += (cvalue - co_offset) * (points - idelay - duration / (0.002 * dhd1.ideltt));
	            if (cvalue > cvmax)
		        cvmax = cvalue;
	            if (cvalue < cvmin)
		        cvmin = cvalue;
		}
	    }
	}
	points++;
	if (points >= data_size + idelay)
	    break;
    }
    points -= idelay;
    if (points != data_size) {
	printf("otto: error in data sequence:\n");
	printf("   points=%d vs expected %d\n", points, data_size);
	printf("      --possible time gap?\n");
    }
    if (dpoints <= 0) {
	printf(" ERROR, dpoints =%d\n", dpoints);
	printf(" (maybe an erroneous variable name)\n");
	return;
    }
    tas_average /= dpoints;	/* doesn't work, accum. mode */
    mean /= dpoints;
    mean += offset;
    trend /= (dpoints / 12. * dpoints * dpoints);
    co_mean /= dpoints;
    co_mean += co_offset;
    co_trend /= (dpoints / 12. * dpoints * dpoints);
    if (ishowcmem) { printf(" cvmax is %f\n", cvmax);}
    				/* for fft case, pad data series as needed */
    if (pmode == 1) {
	while (points < buffer_size) {
	    value = mean + trend * (points - midpoint);
	    *(data_buffer + points) = value;
	    if (ishowcfft) {
		cvalue = co_mean + co_trend * (points - midpoint);
		*(co_data_buffer + points) = cvalue;
	    }
	    points++;
	}
    }
    data_index = 0;
    return;
}

void
plotdata(data_b, mean_local, trend_local, variable_local, proj_local,
	 ylw, yhgh)
    float          *data_b;
    double          mean_local, trend_local;
    char            variable_local[];
    int             proj_local;
    float          *ylw, *yhgh;

{
    float           t, midpoint;
    char            label[46];
    float           ylow = *ylw, yhigh = *yhgh;
    int		    gap = FALSE;

    if ((-1. * ylow) > yhigh)
	yhigh = -ylow;
    if ((-1. * yhigh) < ylow)
	ylow = -1. * yhigh;

    ncar_(&zero, &duration, &ylow, &yhigh, &itwo, &ifive, &itwo, &ifive, &ione, &ione);
    lablx_("TIME :L:(S):U:", 14);
    lably_(variable_local, strlen(variable_local));
    gcolor_(&icolor);

    t = 0.;
    data_index = 1;
    value = *data_b;
    frstpt_(&t, &value);
    while (data_index < buffer_size) {
	t += 0.001 * dhd1.ideltt;
	value = *(data_b + data_index++);
	if (value != MISSING_DATA) {
	    if (gap) {
	        frstpt_(&t, &value);
		gap = FALSE;
	    } else 
	        vector_(&t, &value);
	} else
	    gap = TRUE;
	    
    }

/* rewind and repeat, removing mean and trend */
    gcolor_(&ithree);

    t = 0.;
    data_index = 1;
    gap = FALSE;  
//  midpoint = ((duration - 1.) / (0.002 * dhd1.ideltt));
    midpoint = ((duration) / (0.002 * dhd1.ideltt));
    value = *data_b - (mean_local - trend_local * midpoint);
    frstpt_(&t, &value);
    gap = FALSE;
    while (data_index < buffer_size) {
	t += 0.001 * dhd1.ideltt;
	if ((value = *(data_b + data_index)) != MISSING_DATA) {
	    value -= (mean_local + trend_local * (data_index - midpoint));
	    if (gap) {
		frstpt_(&t, &value);
		gap = FALSE;
	    } else
	        vector_(&t, &value);
	} else
	    gap = TRUE;
        data_index++;
    }
    data_index = 0;
    head_(&date, &start_time, &end_time, &dummy, 1);
    titleline = 1;
    fnote_(&titleline, variable_local, strlen(variable_local));
    if (!strncmp(variable_local, "WIND", 4)) {
	titleline = 3;
	sprintf(label, "wind in %d degree direction", proj_local);
	idline_(&titleline, label, strlen(label));
    }
    gcolor_(&ione);

    /*
     * this strange call forces color reset, needed if plot frames are
     * overlaid by 'med'
     */
    line_(&thousand, &thousand, &thousand, &thousand);
    frame_();
}

void
plotfcst(data_b, mean_local, trend_local, variable_local, proj_local,
	 ylw, yhgh)
    float          *data_b;
    double          mean_local, trend_local;
    char            variable_local[];
    int             proj_local;
    float          *ylw, *yhgh;

{
    char            pyvar[NSZ];
    float           t, midpoint;
    float	    fduration = 1.1 * duration;
    char            label[46];
    float           ylow = *ylw, yhigh = *yhgh;
    float	    size = 10., angle = 0., center = 0.;
    int		    gap = FALSE;
    int		    il;
    int		    iyr, refdays;
    int             izero = 0, ione = 1, itwo = 2, ithree = 3, ifour = 4, 
		    ifive = 5, iseven = 7, ieight = 8, iten = 10;

/* old, from plotdata; suppress: */
#if(0)
    if ((-1. * ylow) > yhigh)
	yhigh = -ylow;
    if ((-1. * yhigh) < ylow)
	ylow = -1. * yhigh;
#endif
    if (variable_local[0] == 'L') {
		/* try to adjust to even values in logarithm */
        il = (10.*yhigh+1.);
	yhigh = il/10.;
        il = 10.*ylow;
	ylow = il/10.;
    } else {
        if (yhigh > 0.) 
            yhigh *= 1.2;
        else
	    yhigh *= 0.8;
        if (ylow > 0.)
            ylow *= 0.8;
        else
	    ylow *= 1.2;
    }
    ncar_(&zero, &fduration, &ylow, &yhigh, &itwo, &ifive, &itwo, &ifive, &ione, &ione);
    lablx_("TIME :L:(S):U:", 14);
    lably_(variable_local, strlen(variable_local));
// draw lines marking months and years (approx):
//	(using 21 market days per month)
    gcolor_(&iten);
    setusv_("LW", &thousand, 2);
    t = duration;
// 	only for durations of less than 10 y:
    if (duration < 21*12*10) {
	iyr = dhd1.idated[0];
	refdays = (iyr+10) * 21 * 12 ;
	refdays += dhd1.idated[2] + 21 * (dhd1.idated[1]-1);
	t -= (float) (refdays % 63);
	refdays -= (refdays % 63);
	do {
	    if ((refdays % 252) == 0) {
	         setusv_("LW", &three_thousand, 2);
	    } else if ((refdays % 126) == 0) {
	         setusv_("LW", &two_thousand, 2);
	    } else {
	         setusv_("LW", &thousand, 2);
	    }
	    frstpt_(&t, &ylow);
	    vector_(&t, &yhigh);
	    t -= 63.;
	    refdays -= 63;
	} while (t > 0.) ;
    }
    gcolor_(&ithree);
    setusv_("LW", &two_thousand, 2);
    frstpt_(&duration, &ylow);
    vector_(&duration, &yhigh);
		// add trend line
//  midpoint = ((duration - 1.) / (0.002 * dhd1.ideltt));
    midpoint = ((duration) / (0.002 * dhd1.ideltt));
    value = mean_local - trend_local * midpoint;
    frstpt_(&zero, &value);
    value = mean_local + trend_local * midpoint * 1.2;
    vector_(&fduration, &value);
		// and 2-sigma band
    setusv_("LW", &thousand, 2);
    LogRMS = (pow(10., LogRMS) - 1.) * 100.;
    value = mean_local - trend_local * midpoint - LogRMS/100.;
    frstpt_(&zero, &value);
    value = mean_local + trend_local * midpoint * 1.2 - LogRMS/100.;
    vector_(&fduration, &value);
    value = mean_local - trend_local * midpoint + LogRMS/100.;
    frstpt_(&zero, &value);
    value = mean_local + trend_local * midpoint * 1.2 + LogRMS/100.;
    vector_(&fduration, &value);
    gcolor_(&icolor);
    setusv_("LW", &three_thousand, 2);

    t = 0.;
    data_index = 1;
    gap = FALSE;  
    value = *data_b + (mean_local - trend_local * midpoint);
    frstpt_(&t, &value);
    gap = FALSE;
    while (data_index < fbuffer_size) {
	t += 0.001 * dhd1.ideltt;
//
//if ((value = *(data_b + data_index)) != MISSING_DATA) {
//              in detrend, missing-value points set to zero
	if ((value = *(data_b + data_index)) != 0.) {
	    value += (mean_local + trend_local * (data_index - midpoint));
	    if (gap) {
		frstpt_(&t, &value);
		gap = FALSE;
	    } else
	        vector_(&t, &value);
	} else
	    gap = TRUE;
        data_index++;
    }
    data_index = 0;
    head_(&date, &start_time, &end_time, "FORECAST", 8);
    titleline = 1;
    setusv_("LW", &thousand, 2);
// add forecast points:
    gcolor_(&iseven);
    t = duration + pday1;
    value = log10((double) value_Day1);
    plchhq_(&t, &value, "*", &size, &angle, &center, 1);
    t = duration + pday2;
    value = log10((double) value_Day2);
    plchhq_(&t, &value, "*", &size, &angle, &center, 1);
    gcolor_(&ione);
    setusv_("LW", &thousand, 2);
#if(STOCKS)
//	now add graph of value of position held:
    if (cvmax < 2500.)
	cvmax = 2500.;
    else if (cvmax < 5000.)
	cvmax = 5000.;
    else if (cvmax < 10000.)
	cvmax = 10000.;
    else if (cvmax < 15000.)
	cvmax = 15000.;
    else if (cvmax < 20000.)
	cvmax = 20000.;
    else if (cvmax < 25000.)
	cvmax = 25000.;
    else
	cvmax = 50000.;
    set_(&ncg.xl1, &ncg.xl2, &ncg.yl1, &ncg.yl2,&zero, &fduration, &zero, &cvmax, &izero);
    axisr_(&izero, &ione, &ione, &ifive);
    gcolor_(&ifour);
    t = duration - 0.001 * dhd1.ideltt;
    data_index = 1;
    gap = FALSE;  
    value = *co_data_buffer;
    CurrentValue = value;
    frstpt_(&t, &value);
    gap = FALSE;
    while (data_index < buffer_size) {
	t -= 0.001 * dhd1.ideltt;
	if ((value = *(co_data_buffer + data_index)) != MISSING_DATA) {
//	    value += (mean_local + trend_local * (data_index - midpoint));
	    if (gap) {
		frstpt_(&t, &value);
		gap = FALSE;
	    } else
	        vector_(&t, &value);
	} else
	    gap = TRUE;
        data_index++;
    }
    setusv_("LW", &two_thousand, 2);
    gcolor_(&iseven);
    t = duration - 0.001 * dhd1.ideltt;
    data_index = 1;
    gap = FALSE;  
    value = *cc_data_buffer;
    CurrentCost = value;
    frstpt_(&t, &value);
    gap = FALSE;
    while (data_index < buffer_size) {
	t -= 0.001 * dhd1.ideltt;
	if ((value = *(cc_data_buffer + data_index)) != MISSING_DATA) {
//	    value += (mean_local + trend_local * (data_index - midpoint));
	    if (gap) {
		frstpt_(&t, &value);
		gap = FALSE;
	    } else
	        vector_(&t, &value);
	} else
	    gap = TRUE;
        data_index++;
    }
    setusv_("LW", &thousand, 2);
    gcolor_(&ione);
    Profit = 100. * (CurrentValue - CurrentCost) / CurrentCost;
    sprintf(labelF1, "Latest Price is %6.2f, Value=%5.0f, Profit=%5.0f\%\%", value_Last, CurrentValue, Profit);
    if(variable_local[0] == 'L') {
	titleline = 3;
		/* estimate annual percentage growth from trend */
	value = (pow(10., 250.*trend_local)-1.)*100.;
	printf("for %s, annual growth (linear trend)=%3.0f\%, RMS is %4.1f\%\n", variable_local+1, value, LogRMS);
	sprintf(label, "annual growth trend = %3.0f\%\%, RMS=%4.0f\%\% ", value, LogRMS);
	idline_(&titleline, label, strlen(label));
	titleline = 2;
	idline_(&titleline, labelF1, strlen(labelF1));
	titleline = 1;
	idline_(&titleline, labelF2, strlen(labelF2));
    }

#if(PYPLT)
    if (pyplot) {
//	now generate pylab plot:
    FILE 	    *templatefile, *outfile;
// set up for python-generated plot
// using 3 files: template, OttoPlot.py, and data. 
    char    Template[120];
    char    PythonPlot[120];
    char           *ch;
    strncpy (pyvar, variable_local, NSZ);
    printf ("vl = %s, pyvar = %s\n", variable_local, pyvar);
    if ((ch = memchr (pyvar, ' ', NSZ)) != NULL) {
	printf (" return is %d\n", ch);
	*ch = '\0';
    }
    strcpy(Template, (char *) getenv("XANADU"));
    strcat(Template, "/src/otto/OttoPythonTemplate.py");
    printf (" Template file is %s.\n", Template);
    templatefile = fopen(Template, "r");
    sprintf (PythonPlot, "%sPlot.py", pyvar+1);
    outfile = fopen(PythonPlot, "w");
// get number of variables (numvars):
    int nvr, ipy, jpy;
    nvr = 4;	// Time, price, value, cost
//	create the data file:
    t = 0.;
    dhd1out.nwords = 8;
    int i;
    for (i = 0; i < fbuffer_size; i++) {
        t += 0.001 * dhd1.ideltt;
        ddataout.ihr = (int) (t+0.001) / 3600;
        ddataout.imin = (int) ((t+0.001)-3600.*ddataout.ihr) / 60;
        ddataout.isec = ((int) (t+0.001)) % 60;
        ddataout.imsec = 0;
        ddataout.values[0] = t;
	if ((value = *(data_b + i)) != 0.) {
	    value += (mean_local + trend_local * (i - midpoint));
	} else {value = MISSING_DATA;}
        ddataout.values[1] = value;
	if (i < buffer_size) {
            ddataout.values[2] = *(co_data_buffer + i);
            ddataout.values[3] = *(cc_data_buffer + i);
	} else {
	    ddataout.values[2] = MISSING_DATA;
            ddataout.values[3] = MISSING_DATA;
	}
        (void) fileout (4);
    }
    (void) fileout (5);	// close the file

// Transfer from template to .py up to break line
    char *rdbuf;
    rdbuf = (char *) malloc(101 * sizeof(char));
    while (fgets(rdbuf, 100, templatefile) != NULL) {
         if (strncmp("XX-----BBBBB-----XX", rdbuf, 19)) {
	    fputs (rdbuf, outfile);
         } else {
            break;
         }
    }
    fprintf (outfile, "numvars = %d\n", nvr);
    fprintf (outfile, "Date = %d\n", date);
//	unlike NCAR-g plot, python plot keeps full resolution
    float forecast_length = fduration * 1.2 / 1.1;
    fprintf (outfile, "DLEN = %d\n", (int) forecast_length * 1000 / dhd1.ideltt);
    fprintf (outfile, "Var = np.empty ([%d, DLEN], 'f')\n", nvr);
    fprintf (outfile, "VarName = [None]*%d\n", nvr);
    while (fgets(rdbuf, 100, templatefile) != NULL) {
         if (strncmp("XX-----BBBBB-----XX", rdbuf, 19)) {
            fputs (rdbuf, outfile);
         } else {
            break;
         }
    }
    fprintf (outfile, "DataFile = open('./DXDATA.%s', 'rb')\n", pyvar+1);
    printf ("ch=%d, vl = %s, pyvar = %s.\n", ch, variable_local, pyvar);
    fprintf (outfile, "VarName[1] = \'%s\'\n", pyvar);
    fprintf (outfile, "%s = np.empty (DLEN, 'f')\n", pyvar);
    pyvar[0] = 'N';
    fprintf (outfile, "VarName[2] = \'%s\'\n", pyvar);
    fprintf (outfile, "%s = np.empty (DLEN, 'f')\n", pyvar);
    pyvar[0] = 'C';
    fprintf (outfile, "VarName[3] = \'%s\'\n", pyvar);
    fprintf (outfile, "%s = np.empty (DLEN, 'f')\n", pyvar);
    fprintf (outfile, "Symbol = VarName[1][1:]\n", pyvar);
    pyvar[0] = 'L';
//    fprintf (outfile, 
    while (fgets(rdbuf, 100, templatefile) != NULL) {
         if (strncmp("XX-----BBBBB-----XX", rdbuf, 19)) {
            fputs (rdbuf, outfile);
         } else {
            break;
         }
    }
    fprintf (outfile, "    Vbuf = struct.unpack (%d*'f', DataFile.read ((numvars)*4))\n", nvr);
    fprintf (outfile, "    Time[i] = Vbuf[0]\n");
    fprintf (outfile, "    %s[i] = Vbuf[1]\n", pyvar);
    pyvar[0] = 'N';
    fprintf (outfile, "    %s[i] = Vbuf[2]\n", pyvar);
    pyvar[0] = 'C';
    fprintf (outfile, "    %s[i] = Vbuf[3]\n", pyvar);
    pyvar[0] = 'L';
    printf ("bug pt 2\n");
    fprintf (outfile, "%s = np.ma.masked_where (%s == -32767., %s)\n", pyvar, pyvar, pyvar);
    printf ("bug pt 3\n");
    pyvar[0] = 'N';
    fprintf (outfile, "%s = np.ma.masked_where (%s == -32767., %s)\n", pyvar, pyvar, pyvar);
    printf ("bug pt 4\n");
    pyvar[0] = 'C';
    fprintf (outfile, "%s = np.ma.masked_where (%s == -32767., %s)\n", pyvar, pyvar, pyvar);
    printf ("bug pt 5\n");
    pyvar[0] = 'L';
    while (fgets(rdbuf, 100, templatefile) != NULL) {
         if (strncmp("XX-----BBBBB-----XX", rdbuf, 19)) {
            fputs (rdbuf, outfile);
         } else {
            break;
         }
    }
    fprintf (outfile, "%s = 10.**%s\n", pyvar, pyvar);
    fprintf (outfile, "LVAR = %s\n", pyvar);
    fprintf (outfile, "Fig = pylab.figure(figsize=(14.,14.))\n");
    fprintf (outfile, "cid = Fig.canvas.mpl_connect ('button_press_event', onclick)\n");
    fprintf (outfile, "Fig.patch.set_facecolor ('#dddddd')\n");
    fprintf (outfile, "PanelR = Fig.add_subplot (111, axisbg='#ccddee')\n");
    fprintf (outfile, "PanelR.yaxis.tick_right()\n");
    fprintf (outfile, "PanelR.yaxis.set_label_position('right')\n");
    fprintf (outfile, "pylab.xlim([0.,%f])\n", forecast_length);
    fprintf (outfile, "pylab.plot(%d-indx, N%s, color='darkgreen')\n", buffer_size, pyvar+1);
    fprintf (outfile, "pylab.plot(%d-indx, C%s, color='orange')\n", buffer_size, pyvar+1);
    fprintf (outfile, "pylab.ylabel('Investment ($)')\n");
    fprintf (outfile, "pylab.plot ([%d,%d], PanelR.get_ylim(), color='purple', lw=2)\n", buffer_size, buffer_size);
// 		// section for display of cursor values
        fprintf (outfile, "def format_coord (x, y):\n");
        fprintf (outfile, "    global DLEN\n");
        fprintf (outfile, "    t = int (x-%d)\n", buffer_size);
        fprintf (outfile, "    if x < DLEN:\n");
        fprintf (outfile, "        y1 = %s[x]\n", pyvar);
	pyvar[0] = 'N';
        fprintf (outfile, "        y2 = %s[%d-x]\n", pyvar, buffer_size);
	pyvar[0] = 'C';
        fprintf (outfile, "        y3 = %s[%d-x]\n", pyvar, buffer_size);
	pyvar[0] = 'L';
        fprintf (outfile, "    else:\n");
        fprintf (outfile, "        y1 = -32767.\n");
        fprintf (outfile, "        y2 = -32767.\n");
        fprintf (outfile, "        y3 = -32767.\n");
        fprintf (outfile, "    return '%%d (%%.2f, %%.0f, %%.0f)'%%(t,y1,y2,y3)\n");
    fprintf (outfile, "PanelL = Fig.add_subplot (111, sharex=PanelR, frameon=False)\n");
    fprintf (outfile, "pylab.plot (indx, %s, label = '%s', color='red', lw=2)\n", pyvar, pyvar+1);
    fprintf (outfile, "pylab.ylabel ('%s price', color='red', fontsize=20)\n", pyvar+1);
    fprintf (outfile, "PanelL.annotate ('latest price %.2f, position %.0f, profit %.0f\%\\nannual trend %3.0f\% RMS %4.0f\%\\nForecasts: %d-day=%.2f; %d-day=%.2f', xycoords='axes fraction', xy=(0.05,0.95), backgroundcolor='yellow')\n", value_Last, CurrentValue, 100.*(CurrentValue-CurrentCost)/CurrentCost, (pow(10., 252.*trend_local)-1.)*100., LogRMS, pday1, value_Day1, pday2, value_Day2);
    fprintf (outfile, "PanelL.annotate ('profit code:\\n gray < 50\%\\n black 50-100\\n orange 100-500\\n yellow 500-1000\\n red 1000-2000\\n white > 2000\%', xycoords='axes fraction', xy=(0.85,0.05), backgroundcolor='yellow')\n");
    fprintf (outfile, "PanelL.set_yscale('log')\n");
    fprintf (outfile, "formatter = pylab.FuncFormatter(log_10_product)\n");
    fprintf (outfile, "PanelL.yaxis.set_major_formatter(formatter)\n");
    fprintf (outfile, "PanelL.yaxis.set_minor_formatter(formatter)\n");
    fprintf (outfile, "if ymin > 0. and ymax > 0.: pylab.ylim ([ymin, ymax])\n");
    fprintf (outfile, "pylab.xlabel('Time [days]')\n");
    fprintf (outfile, "pticks = []\n");
    fprintf (outfile, "pticklabels = []\n");
    t = duration;
// 	only for durations of less than 10 y:
    float lwdth;
    if (duration < 21*12*10) {
	iyr = dhd1.idated[0] + 2000;
	t -= (21./30.) * dhd1.idated[2] + 21 * (dhd1.idated[1]-1);
	t += 2*12*21;
	int pycycle = -1;
	do {
	    pycycle += 1;
	    if (pycycle%4 == 0) {
	         lwdth=2.;
		 if (t < fduration) {
		     fprintf (outfile, "pticks.append (%d)\n", (int) t);
		     fprintf (outfile, "pticklabels.append ('%d')\n", iyr+2- pycycle/4);
		 }
	    } else if (pycycle%4 == 2) {
	         lwdth=1.;
	    } else {
	         lwdth = 0.5;
	    }
	    fprintf (outfile, "pylab.plot ([%d,%d], PanelL.get_ylim(), color='k', lw=%f)\n", (int) t, (int) t, lwdth);
	    t -= 63.;
	} while (t > 0.) ;
    }
    fprintf (outfile, "PanelL.xaxis.set_ticks(pticks)\n");
    fprintf (outfile, "PanelL.xaxis.set_ticklabels(pticklabels)\n");
    fprintf (outfile, "ylim = PanelL.get_ylim()\n");
    fprintf (outfile, "pylab.ylim(ylim)\n");
    fprintf (outfile, "pylab.plot ([0,%d], [%f,%f], color='cyan',lw=1.5)\n", (int) forecast_length, (float) pow(10.,(mean_local-trend_local*midpoint)), (float) pow(10.,(mean_local+trend_local*(forecast_length-midpoint))));
    fprintf (outfile, "pylab.plot ([0,%d], [%f,%f], color='cyan',lw=1.0)\n", (int) forecast_length, (float) pow(10.,(mean_local-trend_local*midpoint - LogRMS/100.)), (float) pow(10.,(mean_local+trend_local*(forecast_length-midpoint) - LogRMS/100.)));
    fprintf (outfile, "pylab.plot ([0,%d], [%f,%f], color='cyan',lw=1.0)\n", (int) forecast_length, (float) pow(10.,(mean_local-trend_local*midpoint + LogRMS/100.)), (float) pow(10.,(mean_local+trend_local*(forecast_length-midpoint) + LogRMS/100.)));
		// read the option-processing code from the template
    while (fgets(rdbuf, 100, templatefile) != NULL) {
         if (strncmp("XX-----BBBBB-----XX", rdbuf, 19)) {
            fputs (rdbuf, outfile);
         } else {
            break;
         }
    }
    free (rdbuf);
    fprintf (outfile, "PanelL.format_coord = format_coord\n");
//    fprintf (outfile, "datacursor (draggable=True, formatter='{label}'.format)\n");
    fprintf (outfile, "pylab.savefig('%sPlot.png')\n", pyvar+1);
    fprintf (outfile, "if show_plot: # pylab.show ()\n");
    printf (" ready to close outfile\n");
    (void) fclose(outfile);
    char syscmd[50];
    sprintf (syscmd, "mv DXDATO DXDATA.%s", pyvar+1);
    (void) system (syscmd);
    if (!strncmp ("GSPC", pyvar+1, 4)) {
        sprintf (syscmd, "python %sPlot.py 0 0 1000. 2000. &", pyvar+1); // flags: plot, download option info
    } else if (!strncmp ("GLD", pyvar+1, 3)) {
        sprintf (syscmd, "python %sPlot.py 0 0 100. 200. &", pyvar+1); // flags: plot, download option info
    } else {
        sprintf (syscmd, "python %sPlot.py 0 0&", pyvar+1); // flags: plot, download option info
    }
    (void) system (syscmd);
#endif
    }  // terminates the python-plot section
    /*
     * this strange call forces color reset, needed if plot frames are
     * overlaid by 'med'
     */
#endif
    line_(&thousand, &thousand, &thousand, &thousand);
    frame_();
}


void
acv(accumulate)
    int             accumulate;	/* 0 for interm., 1 for last */

{

    int            *nr;
    float          *rsave;
    float           dmx;
    float           x, y;
    double          dzero = 0.;
    int             kount = 0, np, i, i2, inda, inds;
    int             ipc, ipcold = 20;
    float           scale_length = 0.;
    float           dx;
    char            label[13];
    static int      acc = 0;
    static int      sum_points = 0;

	maxlags = maxlags_s * sample_frequency;
	plotlags = plotlags_s * sample_frequency;
	smoothpts = smooth_s * sample_frequency;
    if (maxlags >= points) {
	if (accumulate != 0 || acc == 1) {
	    printf(" points must exceed maxlags in accumulate mode; segment skipped\n");
	    fprintf(stderr, " points must exceed maxlags in accumulate mode; segment skipped\n");
	    return;
	}
	max = points - 1;
    } else
	max = maxlags;
    printf(" attempt #9 to allocate %d float words\n", maxlags+1);
    rsave = (float *) GetMemory((size_t) (maxlags+1) * sizeof(float));
    if (acc == 0) {
    printf(" attempt #10 to allocate %d double words\n", maxlags+1);
	r = (double *) GetMemory((size_t) (maxlags+1) * sizeof(double));
    printf(" attempt #11 to allocate %d int words\n", maxlags+1);
	nr = (int *) GetMemory((size_t) (maxlags+1) * sizeof(int));
	acc = 1;
	for (i = 0; i < maxlags; i++) {
	    rsave[i] = 0.;
	    r[i] = dzero;
	    nr[i] = 0;
	}
    }
    data_index = 0;

    while (data_index < data_size) {
	value = *(data_buffer + data_index++);
	inds = kount % maxlags;
	if (value != MISSING_DATA) {
	    rsave[inds] = value;
	    if (maxlags <= kount)
	        i2 = maxlags;
	    else
	        i2 = kount;
	    for (i = 0; i < i2; i++) {
	        inda = inds - i;
	        if (inda < 0)
		    inda += maxlags;
	        r[i] += value * rsave[inda];
		nr[i]++;
	    }
	} else
	    rsave[inds] = 0.; 
	sum_points++;
	if (!(++kount % 100)) {
	    ipc = kount * 40 / points + 20;
	    if (ipc != ipcold) {
		ipcold = ipc;
	    }
	}
    }

    free(rsave);
    if (accumulate)
	return;
    free(data_buffer);
    acc = 0;

    variance = r[0] / sum_points;
    variance = r[0] / nr[0];	// change 110819
    for (i = 1; i < max; i++)
	r[i] /= nr[i] * variance;	// change 110819
//	r[i] /= (r[0]);		/* convert to autocorrelation */
    r[0] = 1.;
    sum_points = 0;
    free(nr);

    if (ishowacv&16) {
        dmx = plotlags_s * tas_average * 0.001;
        printf(" dmx=%f, plotlags_s=%d\n", dmx, plotlags_s);
        ncar_(&zero, &dmx, &mone, &one, &ifive, &itwo, &itwo, &ifive, &ione, &ione);
        lablx_("DISTANCE :L:(KM)", 16);
        lably_("AUTOCORRELATION", 15);
        gcolor_(&icolor);
        x = 0.;
        y = 1.;
        frstpt_(&x, &y);
        dx = tas_average * 0.001 / sample_frequency;
        for (i = 1; i < max; i++) {
	    x = i * dx;
	    y = r[i];
	    scale_length += (y + r[i - 1]) / 2. * dx;
	    vector_(&x, &y);
	    if (i >= plotlags)
	        break;
        }

        sprintf(label, "L=%7.2f km", scale_length);
        printf(" integral scale length = %f\n", scale_length);
        scale_length = 0.;

        if (smooth_s < 1800 && smooth_s > 10) {
	    gcolor_(&ithree);	/* repeat with smoothing */
	    x = 0.;
	    y = 1.;
	    frstpt_(&x, &y);
	    for (i = 1; i < max; i++) {
	        x = i * dx;
	        if (i < smoothpts)
		    r[i] *= (1. - ((float) i) / smoothpts);
	        else
		    r[i] = 0.;
	        y = r[i];
	        scale_length += (y + r[i - 1]) / 2. * dx;
	        vector_(&x, &y);
	    }
	    gcolor_(&icolor);
        }
        printf(" smoothed integral scale length = %f\n", scale_length);
        head_(&date, &start_time, &end_time, &dummy, 1);
        titleline = 1;
        fnote_(&titleline, variable, strlen(variable));
        titleline = 2;
        fnote_(&titleline, "ACV", 3);
        titleline = 2;
        idline_(&titleline, label, 12);
        if (!strncmp(variable, "WIND", 4)) {
	    titleline = 3;
	    sprintf(label, "wind in %d degree direction", proj_degrees);
	    idline_(&titleline, label, strlen(label));
        }
        gcolor_(&ione);

    /*
     * this strange call forces color reset, needed if plot frames are
     * overlaid by 'med'
     */
        line_(&thousand, &thousand, &thousand, &thousand);
        frame_();
        return;
    }
}

void
ftran()
{
    float          *ps, *nu, frequency;
    double          aint, alpha, f2pi;
    float           upl, botl;
    float           flow = .001;
    float           fhigh = 10.;
    float           a, ae, a1, a2;
    double          b;
    int             max2, m, i, j;
    int             label_length;
    int             ipc, ipcold = 60;
    char            label[46];
    float           df, logf;

    max2 = max / 2;
    printf(" attempt #12 to allocate %d float words\n", max2+1);
    ps = (float *) GetMemory((size_t) (max2 + 1) * sizeof(float));
    printf(" attempt #13 to allocate %d float words\n", max2+1);
    nu = (float *) GetMemory((size_t) (max2 + 1) * sizeof(float));

    for (i = 1; i <= max2; i++) {
	frequency = i * (float) sample_frequency / max;
	f2pi = frequency * TWOPI;
	aint = 1.;
	for (j = 1; j < max; j++) {
	    alpha = j / (float) sample_frequency;
	    aint += r[j] * 2. * cos(f2pi * alpha);
	}
	/* aint /= 2.;   -- had this before Aug 1993 */

	/*
	 * should be omitted; 2 above is for -inf--inf integration, 2 below
	 * is for normalization of PSD for freq 0 to inf
	 */

	ps[i] = 2. * aint * variance / (float) sample_frequency;

	nu[i] = frequency;
	ipc = i * 40 / max2 + 60;
	if (ipc != ipcold) {
	    ipcold = ipc;
	}
    }

    if (sample_frequency > 10.) {
	fhigh *= 10.;
    }
    if (ishowacv & 2) {
	i = log10((double) variance) + 3.;
	upl = pow((double) 10., (double) i);
	botl = pow((double) 10., (double) (i - 6));
	if (flow_spec > 0.)
	    flow = flow_spec;
	if (fhigh_spec > 0.)
	    fhigh = fhigh_spec;
	if (plow_spec > 0.)
	    botl = plow_spec;
	if (phigh_spec > 0.)
	    upl = phigh_spec;
	ncar_(&flow, &fhigh, &botl, &upl, &ione, &ione, &ione, &ione, &ifour, &ione);
	lablx_(":PGL:T:RU: :L:(:U:H:L:Z):U:", 27);
	    ilbl = 0;
	    if (!strncmp(variable, "WI", 2)) ilbl = 1;
	    if (!strncmp(variable, "LONW", 4)) ilbl = 1;
	    if (!strncmp(variable, "LATW", 4)) ilbl = 1;
	    if (ilbl)
	        lably_("P(:GL:T:RU:) [m:S1:2 s:S2:-1  ]", 31);
	    else
	        lably_("P(:GL:T:RU:)", 12);
	add_lscale(&flow, &fhigh, &botl, &upl);
	gcolor_(&icolor);
	/* now average together values within smooth_bins */
	df = (log10((double) fhigh) - log10((double) flow)) / smooth_bins;
	logf = log10((double) nu[2]);
	frstpt_(nu + 1, ps + 1);
	i = 1;
	a = 0.;
	frequency = 0.;
	m = 0;
	while (i < max2) {
	    i++;
	    while (log10((double) nu[i]) - logf > df) {
		if (m) {
		    a /= m;
		    frequency /= m;
		    if (a < botl)
			a = botl;
		    vector_(&frequency, &a);
		    a = frequency = 0.;
		    m = 0;
		}
		logf += df;
	    }
	    frequency += nu[i];
	    a += ps[i];
	    m++;
	    if (smooth_bins == MAXSMOOTH) {
		vector_(&frequency, &a);
		frequency = a = 0.;
		m = 0;
	    }
	}
	if (m) {
	    a /= m;
	    frequency /= m;
	    if (a < botl)
		a = botl;
	    vector_(&frequency, &a);
	}
	gcolor_(&ione);
	a = fhigh * 0.9;
	upl *= 0.5;
	cntr = 1.;
	size = 12.;
	angd = 0.;
	(void) sprintf(label, "VARIANCE=%15.7g", variance);
	label_length = 24;
	/* plchhq_(&a, &upl, label, &size, &angd, &cntr, label_length); */
	upl *= 2.;
	head_(&date, &start_time, &end_time, &dummy, 1);
	titleline = 1;
	fnote_(&titleline, variable, strlen(variable));
	idline_(&titleline, label, 24);
	titleline = 2;
	fnote_(&titleline, "ACV", 3);
	sprintf(label, "smoothing interval=%6d seconds", smooth_s);
	idline_(&titleline, label, 33);
	if (!strncmp(variable, "WIND", 4)) {
	    titleline = 3;
	    sprintf(label, "wind in %d degree direction", proj_degrees);
	    idline_(&titleline, label, strlen(label));
	}
	gcolor_(&iseven);
	gsclip_(&ione);
	fhigh /= 10.;
	a = botl * pow((fhigh / flow), 1.666667);
	line_(&flow, &a, &fhigh, &botl);
	fhigh *= 10.;
	gcolor_(&ione);

	/*
	 * this strange call forces color reset, needed if plot frames are
	 * overlaid by 'med'
	 */
	line_(&thousand, &thousand, &thousand, &thousand);
	frame_();
    }
    if (ishowacv & 4) {
	upl = 1.e-10;
	botl = 1.e10;
	for (i = 1; i <= max2; i++) {	/* convert to f * E(f)  */
	    ps[i] *= nu[i];
	    if (ps[i] < botl)
		botl = ps[i];
	    if (ps[i] > upl)
		upl = ps[i];
	}
	if (botl < 1.e-12 && 1.e-12 < upl)
	    botl = 1.e-12;
	a = log10((double) botl);
	if (a < 0)
	    i = a - 1.;
	else
	    i = a;
	botl = pow((double) 10., (double) i);
	if (upl > 1.e6 && 1.e6 > botl)
	    upl = 1.e6;
	a = log10((double) upl);
	if (a > 0)
	    i = a + 1.;
	else
	    i = a;
	upl = pow((double) 10., (double) i);
	if (fplow_spec > 0.)
	    botl = fplow_spec;
	if (fphigh_spec > 0.)
	    upl = fphigh_spec;
	if (flow_spec > 0.) flow = flow_spec;
	if (fhigh_spec > 0.) fhigh = fhigh_spec;
	ncar_(&flow, &fhigh, &botl, &upl, &ione, &ione, &ione, &ione, &ifour, &ione);
	lablx_(":PGL:T:RU: :L:(:U:H:L:Z):U:", 27);
	    if (!strncmp(variable, "WI", 2)) ilbl = 1;
	    if (!strncmp(variable, "LONW", 4)) ilbl = 1;
	    if (!strncmp(variable, "LATW", 4)) ilbl = 1;
	    if (ilbl)
		lably_(":GL:T:RU:P(:GL:T:RU:) [m:S1:2 s:S2:-2  ]", 40);
	    else
		lably_(":GL:T:RU:P(:GL:T:RU:)                   ", 40);
	add_lscale(&flow, &fhigh, &botl, &upl);
	gcolor_(&icolor);
	/* now average together values within smooth_bins */
	df = (log10((double) fhigh) - log10((double) flow)) / smooth_bins;
	logf = log10((double) nu[2]);
	a = *(ps + 1);
	if (a < botl)
	    a = botl;
	frstpt_(nu + 1, &a);
	i = 0;
	a = 0.;
	frequency = 0.;
	m = 0;
	while (i < max2) {
	    i++;
	    while (log10((double) nu[i]) - logf > df) {
		if (m != 0) {
		    a /= m;
		    frequency /= m;
		    if (a < botl)
			a = botl;
		    vector_(&frequency, &a);
		    a = 0.;
		    frequency = 0.;
		    m = 0;
		}
		logf += df;
	    }
	    frequency += nu[i];
	    a += ps[i];
	    m++;
	    if (smooth_bins == MAXSMOOTH) {
		vector_(&frequency, &a);
		frequency = a = 0.;
		m = 0;
	    }
	}
	if (m != 0) {
	    a /= m;
	    frequency /= m;
	    if (a < botl)
		a = botl;
	    vector_(&frequency, &a);
	}
	for (i = 1; i <= max2; i++) {	/* convert back from f * E(f)  */
  	    ps[i] /= nu[i];
	}
	/* add constant-epsilon lines  */
	gcolor_(&iseven);
	setusv_("LW", &thousand, 2);
	dashdb_(&dash_pat2);
	if (spect_variable == 0 || !strncmp (variable, "UX", 2))
	    ae = 0.15;
	else
	    ae = 0.2;
	for (i = -8; i < 0; i++) {
	    if (i == -4) {
		setusv_("LW", &two_thousand, 2);
		dashdb_(&dash_pattern);
	    } else if (i == -3) {
		setusv_("LW", &thousand, 2);
		dashdb_(&dash_pat2);
	    }
	    a = ae * pow((double) 10., (double) i * 0.66667);
	    a *= pow(tas_average, 0.66667);
	    a1 = a * pow((double) flow, -0.66667);
	    if (a1 > upl) {
		a1 = pow((double) 10., (double) i);
		a1 *= tas_average;
		a1 *= pow((double) upl / ae, (double) -1.5);
		if (a1 < fhigh)
		    frstd_(&a1, &upl);
	    } else
		frstd_(&flow, &a1);
	    a2 = a * pow((double) fhigh, -0.66667);
	    if (a2 < botl) {
		a2 = pow((double) 10., (double) i);
		a2 *= tas_average;
		a2 *= pow((double) botl / ae, (double) -1.5);
		if (a2 > flow && a1 < fhigh) {
		    vectd_(&a2, &botl);
		    if (i == -4) {
			a1 = 1.20 * botl;
			a = (log10((double) upl) - log10((double) botl))
			    / (log10((double) fhigh) - log10((double) flow));
			a = -0.66667 / a;
			a = atan((double) a) * 180. / 3.141592;
			size = 8.;
			cntr = 1.;
			plchhq_(&a2, &a1,
			    "1:KGU:V:PRU:10:S2:-4                         ",
				&size, &a, &cntr, 35);
		    }
		}
	    } else if (a1 < fhigh) {
		vectd_(&fhigh, &a2);
		if (i == -4) {
		    a2 *= 1.20;
		    a = (log10((double) upl) - log10((double) botl))
			/ (log10((double) fhigh) - log10((double) flow));
		    a = -0.66667 / a;
		    a = atan((double) a) * 180. / 3.141592;
		    size = 8.;
		    cntr = 1.;
		    plchhq_(&fhigh, &a2,
			    "1:KGU:V:PRU:10:S2:-4                         ",
			    &size, &a, &cntr, 35);
		}
	    }
	}
	setusv_("LW", &two_thousand, 2);
	gcolor_(&ione);
	a = fhigh * 0.9;
	upl *= 0.5;
	cntr = 1.;
	size = 12.;
	angd = 0.;
	(void) sprintf(label, "VARIANCE=%15.7g", variance);
	label_length = 24;
	/* plchhq_(&a, &upl, label, &size, &angd, &cntr, label_length); */
	head_(&date, &start_time, &end_time, &dummy, 1);
	titleline = 1;
	fnote_(&titleline, variable, strlen(variable));
	idline_(&titleline, label, strlen(label));
	titleline = 2;
	fnote_(&titleline, "ACV", 3);
	if (!strncmp(variable, "WIND", 4)) {
	    titleline = 3;
	    sprintf(label, "wind in %d degree direction", proj_degrees);
	    idline_(&titleline, label, strlen(label));
	}
	gcolor_(&ione);

	/*
	 * this strange call forces color reset, needed if plot frames are
	 * overlaid by 'med'
	 */
	line_(&thousand, &thousand, &thousand, &thousand);
	frame_();
    }
    if (ishowacv & 8) {
	i = log10((double) (variance / points)) + 4;
	upl = pow((double) 10., (double) i);
	botl = pow((double) 10., (double) (i - 4));
	upl = 1.e-3;
	botl = 1.e-8;
	if (sample_frequency > 10.)
	    botl /= 10.;
	if (epslow_spec > 0.)
	    botl = epslow_spec;
	if (epshigh_spec > 0.)
	    upl = epshigh_spec;
	ncar_(&flow, &fhigh, &botl, &upl, &ione, &ione, &ione, &ione, &ifour, &ione);
	lablx_(":PGL:T:RU: :L:(:U:H:L:Z):U:", 27);
	lably_("E(:GL:T:RU:) [m:S1:2 s:S2:-3  ]", 31);
	add_lscale(&flow, &fhigh, &botl, &upl);
	gcolor_(&icolor);
	if (spect_variable == 0 || !strncmp (variable, "UX", 2))
	    a = 2.0;
	else
	    a = 1.5;
	a *= pow((double) (TWOPI / tas_average), TWOTHIRDS);
	for (i = 1; i <= max2; i++) {	/* convert to equiv of eddy diss.
					 * rate */
// corrected Oct 2013: should average before taking pow(..., 1.5)
//	    ps[i] = pow((double) ps[i], (double) (1.5)) * a
//		* pow((double) nu[i], 2.5);
	    ps[i] *= a*pow ((double) nu[i], FIVETHIRDS);
	    if (ps[i] < botl)
		ps[i] = botl;
	}
	/* now average together values within smooth_bins */
	df = (log10((double) fhigh) - log10((double) flow)) / smooth_bins;
	logf = log10((double) nu[2]);
	*(ps+1) = pow ((double) (*(ps+1)), 1.5);
	frstpt_(nu + 1, ps + 1);
	i = 0;
	a = 0.;
	frequency = 0.;
	m = 0;
	while (i < max2) {
	    i++;
	    while (log10((double) nu[i]) - logf > df) {
		if (m != 0) {
		    a /= m;
				// mod Oct 2013: now take power
		    a = pow (a, 1.5);
		    frequency /= m;
		    if (a < botl)
			a = botl;
		    vector_(&frequency, &a);
		    a = 0.;
		    frequency = 0.;
		    m = 0;
		}
		logf += df;
	    }
	    frequency += nu[i];
	    a += ps[i];
	    m++;
	    if (smooth_bins == MAXSMOOTH) {
		a = pow (a, 1.5);	// Mod Oct 2013
		vector_(&frequency, &a);
		frequency = a = 0.;
		m = 0;
	    }
	}
	if (m != 0) {
	    a /= m;
	    a = pow (a, 1.5);	// mod Oct 2013
	    frequency /= m;
	    if (a < botl)
		a = botl;
	    vector_(&frequency, &a);
	}
	gcolor_(&ione);
	a = fhigh * 0.9;
	upl *= 0.5;
	cntr = 1.;
	size = 12.;
	angd = 0.;
	(void) sprintf(label, "VARIANCE=%15.7g", variance);
	label_length = 24;
	/* plchhq_(&a, &upl, label, &size, &angd, &cntr, label_length); */
	head_(&date, &start_time, &end_time, &dummy, 1);
	titleline = 1;
	fnote_(&titleline, variable, strlen(variable));
	idline_(&titleline, label, strlen(label));
	titleline = 2;
	fnote_(&titleline, "ACV", 3);
	if (!strncmp(variable, "WIND", 4)) {
	    titleline = 3;
	    sprintf(label, "wind in %d degree direction", proj_degrees);
	    idline_(&titleline, label, strlen(label));
	}
	gcolor_(&ione);

	/*
	 * this strange call forces color reset, needed if plot frames are
	 * overlaid by 'med'
	 */
	line_(&thousand, &thousand, &thousand, &thousand);
	frame_();
    }
    free(ps);
    free(nu);
    free(r);
    return;
}

void
fftplt(dps, sigma, pd_weight, sp_var, var, proj_d, ish)
    double          dps[];
    double          sigma;
    double	    pd_weight;
    int             sp_var;
    char            var[];
    int             proj_d;
    int		    ish;

{

    float          *nu, *ps, frequency;
    float           upl, botl;
    float           ps_plot;
    float           flow = .001;
    float           fhigh = 10.;
    float           a, a1, a2;
    float           ae;
    int             max2, m, i, j;
    int             label_length;
    int             ipc, ipcold = 60;
    int             isgn;
    int             m_eps = 0;
    float           eps, eps_bar = 0., eps2_bar = 0.;
    float           variance_20km = 0.;
    char            label[50];
    float           df, logf;
    double	    cf;

    max2 = segment_length / 2;
    nu = vector(1, max2);	/* allocate storage */
    ps = vector(1, max2);
    variance = 0.;
    nu[1] = 0.;
    if (sp_var == 0 || !strncmp (variable, "UX", 2))
	ae = 2.0;
    else
	ae = 1.5;
    ae *= pow((double)(TWOPI / tas_average), TWOTHIRDS);
    /* needed below: */
    nu[max2] = max2 * ((float) sample_frequency) / segment_length;
    cf = ((double) segment_length) / (double) sample_frequency;
    for (i = 1; i <= max2; i++) {
	nu[i] = i / cf;
//	dps[i] *= pd_weight;	
//  this error was corrected 110927: changing dps changed fps and so
//  affected subsequent calculation of coherence
	variance += dps[i] * pd_weight;
	if (nu[i] > tas_average / 20000.)	/* special KOFSE, temporary */
//	    variance_20km += dps[i];
	    variance_20km += dps[i] * pd_weight;
//	ps[i] = (float) (dps[i] * cf);
	ps[i] = (float) (dps[i] * cf * pd_weight);

	/*
	 * calculate average eddy dissipation rate, last decade of
	 * frequencies
	 */
	if (nu[i] >= nu[max2] / 10.) {
//	    eps = pow((double) ps[i], (double) (1.5)) * ae
//		* pow((double) nu[i], 2.5);
	    eps = ae * ps[i] * pow ((double) nu[i], FIVETHIRDS);  // corr Oct 2013
	    eps_bar += eps;
	    eps2_bar += eps * eps;
	    m_eps++;
	}
	ipc = i * 10 / max2 + 60;
	if (ipc != ipcold) {
	    ipcold = ipc;
	}
    }
    if (m_eps > 0) {
	eps_bar /= m_eps;
	eps2_bar /= m_eps;
	eps_bar = pow ((double) eps_bar, 1.5); // corr Oct 2013
	eps2_bar = (float) sqrt((double) (eps2_bar - eps_bar * eps_bar) / m_eps);
	eps2_bar = pow ((double) eps2_bar, 1.5);  // corr Oct 2013
	printf(" FFT average eddy dissipation rate=%.3e +/- %.3e\n",
	       eps_bar, eps2_bar);
    }
    if (sample_frequency > 10.) {
	fhigh *= 10.;
    }
    if (ish & 2) {	/* plot P(f) */
	i = log10((double) variance) + 3.;
	upl = pow((double) 10., (double) i);
	botl = pow((double) 10., (double) (i - 6));
	if (upl > 1.e20 || botl < 1.e-20) {
	    printf(" plot outside expected limits, low=%f, high=%f\n",
		   botl, upl);
	    return;
	}
	if (flow_spec > 0.)
	    flow = flow_spec;
	if (fhigh_spec > 0.)
	    fhigh = fhigh_spec;
	if (plow_spec > 0.)
	    botl = plow_spec;
	if (phigh_spec > 0.)
	    upl = phigh_spec;
		// write the spectrum to a binary file
//	create the spectrum file:
        float t = 0.;
        dhd1out.nwords = 7;
	for (i = 1; i <= max2; i++) {
            t += 0.001 * dhd1.ideltt;
            ddataout.ihr = (int) (t+0.001) / 3600;
            ddataout.imin = (int) ((t+0.001)-3600.*ddataout.ihr) / 60;
            ddataout.isec = ((int) (t+0.001)) % 60;
            ddataout.imsec = 0;
            ddataout.values[0] = t;
            ddataout.values[1] = nu[i];
            ddataout.values[2] = ps[i];
            (void) fileout (4);
        }
        (void) fileout (5);	// close the file
        char syscmd[50];	// and rename it to avoid overlap
        sprintf (syscmd, "mv DXDATO DXDATP");
        (void) system (syscmd);
//	now generate pylab plot:
        FILE 	    *templatefile, *outfile;
// set up for python-generated plot
// using 3 files: template, OttoFFT.py, and data. 
        char    Template[120];
        char    PythonPlot[120];
        char           *ch;
        char            pyvar[NSZ];
        strncpy (pyvar, var, NSZ);
        if ((ch = memchr (pyvar, ' ', NSZ)) != NULL) {
	    printf (" return is %d\n", ch);
	    *ch = '\0';
        }
        strcpy(Template, (char *) getenv("XANADU"));
        strcat(Template, "/src/otto/OttoFFTTemplate.py");
        printf (" Template file is %s.\n", Template);
        templatefile = fopen(Template, "r");
//      sprintf (PythonPlot, "%sPlot.py", pyvar+1);
//      outfile = fopen(PythonPlot, "w");
        outfile = fopen ("FFTPlot.py", "w");
        char *rdbuf;
        rdbuf = (char *) malloc(101 * sizeof(char));
        while (fgets(rdbuf, 100, templatefile) != NULL) {
            if (strncmp("XX-----BBBBB-----XX", rdbuf, 19)) {
	        fputs (rdbuf, outfile);
            } else { break; }
        }
        fprintf (outfile, "numvars = 3\n");
        fprintf (outfile, "Date = %d\n", date);
        fprintf (outfile, "DLEN = %d\n", max2);
        fprintf (outfile, "Var = np.empty ([3, DLEN], 'f')\n");
        fprintf (outfile, "VarName = [None]*3\n");
        while (fgets(rdbuf, 100, templatefile) != NULL) {
            if (strncmp("XX-----BBBBB-----XX", rdbuf, 19)) {
                fputs (rdbuf, outfile);
            } else {
                break;
            }
        }
        fprintf (outfile, "DataFile = open('./DXDATP', 'rb')\n");
        fprintf (outfile, "VarName[1] = \'nu\'\n");
        fprintf (outfile, "nu = np.empty (DLEN, 'f')\n");
        fprintf (outfile, "VarName[2] = \'ps\'\n");
        fprintf (outfile, "ps = np.empty (DLEN, 'f')\n");
        while (fgets(rdbuf, 100, templatefile) != NULL) {
            if (strncmp("XX-----BBBBB-----XX", rdbuf, 19)) {
                fputs (rdbuf, outfile);
            } else {
                break;
            }
        }
        fprintf (outfile, "    Vbuf = struct.unpack (3*'f', DataFile.read ((numvars)*4))\n");
        fprintf (outfile, "    Time[i] = Vbuf[0]\n");
        fprintf (outfile, "    nu[i] = Vbuf[1]\n");
        fprintf (outfile, "    ps[i] = Vbuf[2]\n");
        fprintf (outfile, "nu = np.ma.masked_where (nu == -32767., nu)\n");
        fprintf (outfile, "ps = np.ma.masked_where (ps == -32767., ps)\n");
        while (fgets(rdbuf, 100, templatefile) != NULL) {
            if (strncmp("XX-----BBBBB-----XX", rdbuf, 19)) {
                fputs (rdbuf, outfile);
            } else {
                break;
            }
        }
        fprintf (outfile, "Fig = pylab.figure(figsize=(6.5,6))\n");
//      fprintf (outfile, "cid = Fig.canvas.mpl_connect ('button_press_event', onclick)\n");
        fprintf (outfile, "Fig.patch.set_facecolor ('#dddddd')\n");
        fprintf (outfile, "PanelL = Fig.add_subplot (111, axisbg='#ccddee')\n");
// 		// section for display of cursor values
        fprintf (outfile, "def format_coord2 (al, y):\n");
        fprintf (outfile, "    global DLEN, tas_average, avenu, aveps\n");
        fprintf (outfile, "    x = tas_average / al\n");
        fprintf (outfile, "    if x >= avenu[0] and x <= avenu[-1]:\n");
        fprintf (outfile, "        y1 = x\n");
        fprintf (outfile, "        y2 = al\n");
        fprintf (outfile, "        ix = 0\n");
        fprintf (outfile, "        while avenu[ix] < x and ix < len (avenu): ix += 1\n");
        fprintf (outfile, "        y3 = aveps[ix]\n");
        fprintf (outfile, "    else:\n");
        fprintf (outfile, "        y1 = -32767.\n");
        fprintf (outfile, "        y2 = -32767.\n");
        fprintf (outfile, "        y3 = -32767.\n");
        fprintf (outfile, "    return '(%%.2f, %%.2f, %%.2e)'%%(y1,y2,y3)\n");
        fprintf (outfile, "if show_ps:\n");
        fprintf (outfile, "    pylab.plot (nu, ps, label = 'PSD', color='red', lw=0.8)\n");
        fprintf (outfile, "    pylab.ylabel (r'%s:  $\\ $P(f)', color='red', fontsize=20)\n", pyvar);
        fprintf (outfile, "tas_average=%f\n", tas_average);
        fprintf (outfile, "lateral = %d\n", sp_var);
	if (flow_spec > 0.) {
            fprintf (outfile, "flow=%f\n", flow_spec);
	} else {
            fprintf (outfile, "flow=%f\n", flow);
	}
        if (fhigh_spec > 0.) {
            fprintf (outfile, "fhigh=%f\n", fhigh_spec);
	} else {
            fprintf (outfile, "fhigh=%f\n", fhigh);
	}
	fprintf (outfile, "fplow=%f\n", botl);
	fprintf (outfile, "fphigh=%f\n", upl);
        fprintf (outfile, "smooth_bins=%d\n", smooth_bins);
        fprintf (outfile, "variance = %.1f\n", variance);
        fprintf (outfile, "EDR = %.5f\n", eps_bar);
        while (fgets(rdbuf, 100, templatefile) != NULL) {
            if (strncmp("XX-----BBBBB-----XX", rdbuf, 19)) {
                fputs (rdbuf, outfile);
            } else {
                break;
            }
        }
        fprintf (outfile, "PanelL.annotate ('Method: FFT; Variance=%.2f m$^2$s$^{-2}$\\npts/seg: %d, bins for smoothing=%d\\nestimated eddy dissipation rate: %.1e m$^2$s$^{-3}$', xycoords='axes fraction', xy=(0.02, 0.05), backgroundcolor='lightyellow')\n", variance, segment_length, smooth_bins, eps_bar);
        fprintf (outfile, "PanelL.annotate ('%d %d--%d', xycoords='axes fraction', xy=(0.5,0.937), backgroundcolor='lightyellow')\n", date, start_time, end_time);
//      fprintf (outfile, "formatter = pylab.FuncFormatter(log_10_product)\n");
//      fprintf (outfile, "PanelL.yaxis.set_major_formatter(formatter)\n");
//      fprintf (outfile, "PanelL.xaxis.set_major_formatter(formatter)\n");
//      fprintf (outfile, "PanelL.yaxis.set_minor_formatter(formatter)\n");
//      fprintf (outfile, "if ymin > 0. and ymax > 0.: pylab.ylim ([ymin, ymax])\n");
//      fprintf (outfile, "ylim = PanelL.get_ylim()\n");
//      fprintf (outfile, "pylab.ylim(ylim)\n");
        free (rdbuf);
        fprintf (outfile, "pylab.grid()\n");
//        fprintf (outfile, "datacursor (draggable=True, formatter='{label}'.format)\n");
        fprintf (outfile, "flim = PanelL.get_xlim()\n");
        fprintf (outfile, "PanelT = PanelL.twiny()\n");
        fprintf (outfile, "PanelT.set_xlabel('wavelength [m]', color='black')\n");
        fprintf (outfile, "pylab.xlim ([tas_average / flim[0], tas_average / flim[1]])\n");
        fprintf (outfile, "PanelT.set_xscale ('log')\n");
        fprintf (outfile, "for tl in PanelT.get_xticklabels():\n");
        fprintf (outfile, "    tl.set_color('black')\n");
        fprintf (outfile, "PanelT.format_coord = format_coord2\n");

        fprintf (outfile, "pylab.savefig('FFTPlot.png', facecolor='#dddddd')\n");
        fprintf (outfile, "# pylab.show ()\n");
        printf (" ready to close outfile\n");
        (void) fclose(outfile);
        sprintf (syscmd, "python FFTPlot.py &"); 
        (void) system (syscmd);
	ncar_(&flow, &fhigh, &botl, &upl, &ione, &ione, &ione, &ione, &ifour, &ione);
	lablx_(":PGL:T:RU: :L:(:U:H:L:Z):U:", 27);
	    ilbl = 0;
	    if (!strncmp(var, "WI", 2)) ilbl = 1;
	    if (!strncmp(var, "LONW", 4)) ilbl = 1;
	    if (!strncmp(var, "LATW", 4)) ilbl = 1;
	    if (ilbl)
	        lably_("P(:GL:T:RU:) [m:S1:2 s:S2:-1  ]", 31);
	    else
	        lably_("P(:GL:T:RU:)", 12);
	add_lscale(&flow, &fhigh, &botl, &upl);
	df = (log10((double) fhigh) - log10((double) flow)) / smooth_bins;
	for (isgn = -1; isgn <= 1; isgn++) {
	    if (isgn && !show_errors)
		continue;
	    if (isgn) {
		gcolor_(&ithree);
		setusv_("LW", &thousand, 2);
	    } else {
		gcolor_(&icolor);
		setusv_("LW", &two_thousand, 2);
	    }
	    /* now average together values within smooth_bins */
	    logf = log10((double) nu[2]);
	    a = ps[1] * (1. + isgn / sqrt(sigma));
	    if (a < botl)
		a = botl;
	    frstpt_(nu + 1, &a);
#if(CORR)
		// write out values for follow-up analyses
	    if(!isgn) {printf(" %13.5f %13.5f\n", *(nu+1),a);}
#endif 
	    i = 1;
	    a = 0.;
	    frequency = 0.;
	    m = 0;
	    while (i++ < max2) {
		while (log10((double) nu[i]) - logf > df) {
		    if (m) {
			a *= (1. + isgn / sqrt(sigma * m)) / m;
			frequency /= m;
#if(CORR)
		// write out values for follow-up analyses
	    if(!isgn) {printf(" %13.5f %13.5f\n", frequency, a);}
#endif 
			if (a < botl)
			    a = botl;
			vector_(&frequency, &a);
			a = frequency = 0.;
			m = 0;
		    }
		    logf += df;
		}
		frequency += nu[i];
		a += ps[i];
		m++;
		if (smooth_bins == MAXSMOOTH) {
		    a *= (1. + isgn / sqrt(sigma));
#if(CORR)
		// write out values for follow-up analyses
	    if(!isgn) {printf(" %13.5f %13.5f\n", frequency, a);}
#endif 
		    vector_(&frequency, &a);
		    frequency = a = 0.;
		    m = 0;
		}
	    }
	    if (m) {
		a *= (1. + isgn / sqrt(sigma * m)) / m;
		frequency /= m;
#if(CORR)
		// write out values for follow-up analyses
	    if(!isgn) {printf(" %13.5f %13.5f\n", frequency, a);}
#endif 
		if (a < botl)
		    a = botl;
		vector_(&frequency, &a);
	    }
	}
	setusv_("LW", &two_thousand, 2);
	gcolor_(&ione);
	a = fhigh * 0.9;
	upl *= 0.5;
	cntr = 1.;
	size = 12.;
	angd = 0.;
	(void) sprintf(label, "VARIANCE=%15.7g", variance);
	label_length = 24;
	/* plchhq_(&a, &upl, label, &size, &angd, &cntr, label_length); */
	upl *= 2.;
	head_(&date, &start_time, &end_time, &dummy, 1);
	titleline = 1;
	fnote_(&titleline, var, strlen(var));
	idline_(&titleline, label, 24);
	titleline = 2;
	fnote_(&titleline, "FFT", 3);
	sprintf(label, "pts/seg=%6d, window=%2d, smooth bins=%5d",
		segment_length, iwindow + 1, smooth_bins);
	idline_(&titleline, label, 44);
	if (!strncmp(var, "WIND", 4)) {
	    titleline = 3;
	    sprintf(label, "wind in %d degree direction", proj_d);
	    idline_(&titleline, label, strlen(label));
	}
	gcolor_(&iseven);
	gsclip_(&ione);
	fhigh /= 10.;
	a = botl * pow((fhigh / flow), 1.666667);
	line_(&flow, &a, &fhigh, &botl);
	fhigh *= 10.;
	gcolor_(&ione);

	/*
	 * this strange call forces color reset, needed if plot frames are
	 * overlaid by 'med'
	 */
	line_(&thousand, &thousand, &thousand, &thousand);
	frame_();
    }
    if (ish & 4) {	/* plot f P(f) */
	upl = 1.e-10;
	botl = 1.e10;
	for (i = 1; i <= max2; i++) {	/* convert to f * E(f)  */
	    ps[i] *= nu[i];
	    if (ps[i] < botl)
		botl = ps[i];
	    if (ps[i] > upl)
		upl = ps[i];
	}
	if (botl < 1.e-12 && 1.e-12 < upl)
	    botl = 1.e-12;
	a = log10((double) botl);
	if (a < 0)
	    i = a - 1.;
	else
	    i = a;
	botl = pow((double) 10., (double) i);
	if (upl > 1.e6 && 1.e6 > botl)
	    upl = 1.e6;
	a = log10((double) upl);
	if (a > 0)
	    i = a + 1.;
	else
	    i = a;
	upl = pow((double) 10., (double) i);
	if (fplow_spec > 0.)
	    botl = fplow_spec;
	if (fphigh_spec > 0.)
	    upl = fphigh_spec;
	if (flow_spec > 0.) flow = flow_spec;
	if (fhigh_spec > 0.) fhigh = fhigh_spec;
		// write the spectrum to a binary file
//	create the spectrum file:
        float t = 0.;
        dhd1out.nwords = 7;
	for (i = 1; i <= max2; i++) {
            t += 0.001 * dhd1.ideltt;
            ddataout.ihr = (int) (t+0.001) / 3600;
            ddataout.imin = (int) ((t+0.001)-3600.*ddataout.ihr) / 60;
            ddataout.isec = ((int) (t+0.001)) % 60;
            ddataout.imsec = 0;
            ddataout.values[0] = t;
            ddataout.values[1] = nu[i];
            ddataout.values[2] = ps[i];
            (void) fileout (4);
        }
        (void) fileout (5);	// close the file
        char syscmd[50];	// and rename it to avoid overlap
        sprintf (syscmd, "mv DXDATO DXDATS");
        (void) system (syscmd);
//	now generate pylab plot:
        FILE 	    *templatefile, *outfile;
// set up for python-generated plot
// using 3 files: template, OttoFFT.py, and data. 
        char    Template[120];
        char    PythonPlot[120];
        char           *ch;
        char            pyvar[NSZ];
        strncpy (pyvar, var, NSZ);
        if ((ch = memchr (pyvar, ' ', NSZ)) != NULL) {
	    printf (" return is %d\n", ch);
	    *ch = '\0';
        }
        strcpy(Template, (char *) getenv("XANADU"));
        strcat(Template, "/src/otto/OttoFFTTemplate.py");
        printf (" Template file is %s.\n", Template);
        templatefile = fopen(Template, "r");
//      sprintf (PythonPlot, "%sPlot.py", pyvar+1);
//      outfile = fopen(PythonPlot, "w");
        outfile = fopen ("FFTPlot.py", "w");
        char *rdbuf;
        rdbuf = (char *) malloc(101 * sizeof(char));
        while (fgets(rdbuf, 100, templatefile) != NULL) {
            if (strncmp("XX-----BBBBB-----XX", rdbuf, 19)) {
	        fputs (rdbuf, outfile);
            } else { break; }
        }
        fprintf (outfile, "numvars = 3\n");
        fprintf (outfile, "Date = %d\n", date);
        fprintf (outfile, "DLEN = %d\n", max2);
        fprintf (outfile, "Var = np.empty ([3, DLEN], 'f')\n");
        fprintf (outfile, "VarName = [None]*3\n");
        while (fgets(rdbuf, 100, templatefile) != NULL) {
            if (strncmp("XX-----BBBBB-----XX", rdbuf, 19)) {
                fputs (rdbuf, outfile);
            } else {
                break;
            }
        }
        fprintf (outfile, "DataFile = open('./DXDATS', 'rb')\n");
        fprintf (outfile, "VarName[1] = \'nu\'\n");
        fprintf (outfile, "nu = np.empty (DLEN, 'f')\n");
        fprintf (outfile, "VarName[2] = \'ps\'\n");
        fprintf (outfile, "ps = np.empty (DLEN, 'f')\n");
        while (fgets(rdbuf, 100, templatefile) != NULL) {
            if (strncmp("XX-----BBBBB-----XX", rdbuf, 19)) {
                fputs (rdbuf, outfile);
            } else {
                break;
            }
        }
        fprintf (outfile, "    Vbuf = struct.unpack (3*'f', DataFile.read ((numvars)*4))\n");
        fprintf (outfile, "    Time[i] = Vbuf[0]\n");
        fprintf (outfile, "    nu[i] = Vbuf[1]\n");
        fprintf (outfile, "    ps[i] = Vbuf[2]\n");
        fprintf (outfile, "nu = np.ma.masked_where (nu == -32767., nu)\n");
        fprintf (outfile, "ps = np.ma.masked_where (ps == -32767., ps)\n");
        while (fgets(rdbuf, 100, templatefile) != NULL) {
            if (strncmp("XX-----BBBBB-----XX", rdbuf, 19)) {
                fputs (rdbuf, outfile);
            } else {
                break;
            }
        }
        fprintf (outfile, "Fig = pylab.figure(figsize=(6.5,6))\n");
//      fprintf (outfile, "cid = Fig.canvas.mpl_connect ('button_press_event', onclick)\n");
        fprintf (outfile, "Fig.patch.set_facecolor ('#dddddd')\n");
        fprintf (outfile, "PanelL = Fig.add_subplot (111, axisbg='#ccddee')\n");
// 		// section for display of cursor values
        fprintf (outfile, "def format_coord2 (al, y):\n");
        fprintf (outfile, "    global DLEN, tas_average, avenu, aveps\n");
        fprintf (outfile, "    x = tas_average / al\n");
        fprintf (outfile, "    if x >= avenu[0] and x <= avenu[-1]:\n");
        fprintf (outfile, "        y1 = x\n");
        fprintf (outfile, "        y2 = al\n");
        fprintf (outfile, "        ix = 0\n");
        fprintf (outfile, "        while avenu[ix] < x and ix < len (avenu): ix += 1\n");
        fprintf (outfile, "        y3 = aveps[ix]\n");
        fprintf (outfile, "    else:\n");
        fprintf (outfile, "        y1 = -32767.\n");
        fprintf (outfile, "        y2 = -32767.\n");
        fprintf (outfile, "        y3 = -32767.\n");
        fprintf (outfile, "    return '(%%.2f, %%.2f, %%.2e)'%%(y1,y2,y3)\n");
        fprintf (outfile, "if show_ps:\n");
        fprintf (outfile, "    pylab.plot (nu, ps, label = 'PSD', color='red', lw=0.8)\n");
        fprintf (outfile, "    pylab.ylabel (r'%s:  f$\\ $P(f) = $-\\lambda P\\ (\\lambda)$', color='red', fontsize=20)\n", pyvar);
        fprintf (outfile, "tas_average=%f\n", tas_average);
        fprintf (outfile, "lateral = %d\n", sp_var);
	if (flow_spec > 0.) {
            fprintf (outfile, "flow=%f\n", flow_spec);
	} else {
            fprintf (outfile, "flow=%f\n", flow);
	}
        if (fhigh_spec > 0.) {
            fprintf (outfile, "fhigh=%f\n", fhigh_spec);
	} else {
            fprintf (outfile, "fhigh=%f\n", fhigh);
	}
	fprintf (outfile, "fplow=%f\n", fplow_spec);
	fprintf (outfile, "fphigh=%f\n", fphigh_spec);
        fprintf (outfile, "smooth_bins=%d\n", smooth_bins);
        fprintf (outfile, "variance = %.1f\n", variance);
        fprintf (outfile, "EDR = %.5f\n", eps_bar);
        while (fgets(rdbuf, 100, templatefile) != NULL) {
            if (strncmp("XX-----BBBBB-----XX", rdbuf, 19)) {
                fputs (rdbuf, outfile);
            } else {
                break;
            }
        }
        fprintf (outfile, "PanelL.annotate ('Method: FFT; Variance=%.2f m$^2$s$^{-2}$\\npts/seg: %d, bins for smoothing=%d\\nestimated eddy dissipation rate: %.1e m$^2$s$^{-3}$', xycoords='axes fraction', xy=(0.02, 0.05), backgroundcolor='lightyellow')\n", variance, segment_length, smooth_bins, eps_bar);
        fprintf (outfile, "PanelL.annotate ('%d %d--%d', xycoords='axes fraction', xy=(0.5,0.937), backgroundcolor='lightyellow')\n", date, start_time, end_time);
//      fprintf (outfile, "formatter = pylab.FuncFormatter(log_10_product)\n");
//      fprintf (outfile, "PanelL.yaxis.set_major_formatter(formatter)\n");
//      fprintf (outfile, "PanelL.xaxis.set_major_formatter(formatter)\n");
//      fprintf (outfile, "PanelL.yaxis.set_minor_formatter(formatter)\n");
//      fprintf (outfile, "if ymin > 0. and ymax > 0.: pylab.ylim ([ymin, ymax])\n");
//      fprintf (outfile, "ylim = PanelL.get_ylim()\n");
//      fprintf (outfile, "pylab.ylim(ylim)\n");
        free (rdbuf);
        fprintf (outfile, "pylab.grid()\n");
//        fprintf (outfile, "datacursor (draggable=True, formatter='{label}'.format)\n");
        fprintf (outfile, "flim = PanelL.get_xlim()\n");
        fprintf (outfile, "PanelT = PanelL.twiny()\n");
        fprintf (outfile, "PanelT.set_xlabel('wavelength [m]', color='black')\n");
        fprintf (outfile, "pylab.xlim ([tas_average / flim[0], tas_average / flim[1]])\n");
        fprintf (outfile, "PanelT.set_xscale ('log')\n");
        fprintf (outfile, "for tl in PanelT.get_xticklabels():\n");
        fprintf (outfile, "    tl.set_color('black')\n");
        fprintf (outfile, "PanelT.format_coord = format_coord2\n");

        fprintf (outfile, "pylab.savefig('FFTPlot.png', facecolor='#dddddd')\n");
        fprintf (outfile, "# pylab.show ()\n");
        printf (" ready to close outfile\n");
        (void) fclose(outfile);
        sprintf (syscmd, "python FFTPlot.py &"); 
        (void) system (syscmd);
	
	ncar_(&flow, &fhigh, &botl, &upl, &ione, &ione, &ione, &ione, &ifour, &ione);
	lablx_(":PGL:T:RU: :L:(:U:H:L:Z):U:", 27);
	    ilbl = 0;
	    if (!strncmp(var, "WI", 2)) ilbl = 1;
	    if (!strncmp(var, "LONW", 4)) ilbl = 1;
	    if (!strncmp(var, "LATW", 4)) ilbl = 1;
	    if (ilbl)
		lably_(":GL:T:RU:P(:GL:T:RU:) [m:S1:2 s:S2:-2  ]", 40);
	    else
		lably_(":GL:T:RU:P(:GL:T:RU:)                   ", 40);
	add_lscale(&flow, &fhigh, &botl, &upl);
	gcolor_(&icolor);
	/* now average together values within smooth_bins */
	df = (log10((double) fhigh) - log10((double) flow)) / smooth_bins;
	for (isgn = -1; isgn <= 1; isgn++) {
	    if (isgn && !show_errors)
		continue;
	    if (isgn) {
		gcolor_(&ithree);
		setusv_("LW", &thousand, 2);
	    } else {
		gcolor_(&icolor);
		setusv_("LW", &two_thousand, 2);
	    }
	    logf = log10((double) nu[2]);
	    a = ps[1] * (1. + isgn / sqrt(sigma));
	    if (a < botl)
		a = botl;
	    frstpt_(nu + 1, &a);
	    i = 1;
	    a = 0.;
	    frequency = 0.;
	    m = 0;
	    while (i++ < max2) {
		while (log10((double) nu[i]) - logf > df) {
		    if (m != 0) {
			frequency /= m;
			a *= (1. + isgn / sqrt(sigma * m)) / m;
			if (a < botl)
			    a = botl;
			vector_(&frequency, &a);
			a = 0.;
			frequency = 0.;
			m = 0;
		    }
		    logf += df;
		}
		frequency += nu[i];
		a += ps[i];
		m++;
		if (smooth_bins == MAXSMOOTH) {
		    a *= (1. + isgn / sqrt(sigma));
		    vector_(&frequency, &a);
		    frequency = a = 0.;
		    m = 0;
		}
	    }
	    if (m != 0) {
		a *= (1. + isgn / sqrt(sigma * m)) / m;
		frequency /= m;
		if (a < botl)
		    a = botl;
		vector_(&frequency, &a);
	    }
	}
	for (i = 1; i <= max2; i++) {	/* convert back from f * E(f)  */
	    ps[i] /= nu[i];
	}
	/* add constant-epsilon and constant-noise lines  */
	gcolor_(&iseven);
	setusv_("LW", &thousand, 2);
	dashdb_(&dash_pat2);
	if (sp_var == 0 || !strncmp (var, "UX", 2) || !strncmp (var, "TASX", 4))
	    ae = 0.15;
	else
	    ae = 0.2;
        printf (" ae for var=%s and sp_var=%d is %f\n", var, sp_var, ae);
	for (i = -8; i < 0; i++) {
	    if (i == -4) {
		setusv_("LW", &two_thousand, 2);
		dashdb_(&dash_pattern);
	    } else if (i == -3) {
		setusv_("LW", &thousand, 2);
		dashdb_(&dash_pat2);
	    }
	    a = ae * pow((double) 10., (double) i * 0.66667);
	    a *= pow(tas_average, 0.66667);
	    a1 = a * pow((double) flow, -0.66667);
	    if (a1 > upl) {
		a1 = pow((double) 10., (double) i);
		a1 *= tas_average;
		a1 *= pow((double) upl / ae, (double) -1.5);
		if (a1 < fhigh)
		    frstd_(&a1, &upl);
	    } else
		frstd_(&flow, &a1);
	    a2 = a * pow((double) fhigh, -0.66667);
	    if (a2 < botl) {
		a2 = pow((double) 10., (double) i);
		a2 *= tas_average;
		a2 *= pow((double) botl / ae, (double) -1.5);
		if (a2 > flow && a1 < fhigh) {
		    vectd_(&a2, &botl);
		    if (i == -4) {
			a1 = 1.20 * botl;
			a = (log10((double) upl) - log10((double) botl))
			    / (log10((double) fhigh) - log10((double) flow));
			a = -0.66667 / a;
			a = atan((double) a) * 180. / 3.141592;
			size = 8.;
			cntr = 1.;
			plchhq_(&a2, &a1,
			    "1:KGU:V:PRU:10:S2:-4                         ",
				&size, &a, &cntr, 35);
		    }
		}
	    } else if (a1 < fhigh) {
		vectd_(&fhigh, &a2);
		if (i == -4) {
		    a2 *= 1.20;
		    a = (log10((double) upl) - log10((double) botl))
			/ (log10((double) fhigh) - log10((double) flow));
		    a = -0.66667 / a;
		    a = atan((double) a) * 180. / 3.141592;
		    size = 8.;
		    cntr = 1.;
		    plchhq_(&fhigh, &a2,
			    "1:KGU:V:PRU:10:S2:-4                         ",
			    &size, &a, &cntr, 35);
		}
	    }
	}
	setusv_("LW", &thousand, 2);
	dashdb_(&dash_pat2);
	for (i = -8; i <= 2; i++) {
	    if (i == -4) {
		setusv_("LW", &two_thousand, 2);
		dashdb_(&dash_pattern);
	    } else if (i == -3) {
		setusv_("LW", &thousand, 2);
		dashdb_(&dash_pat2);
	    }
	    a = pow((double) 10., (double) i);
	    a1 = a * flow;
	    if (a1 < botl) {
		a1 = botl / a;
		if (a1 > flow) {
		    frstd_(&a1, &botl);
		}
	    } else
		frstd_(&flow, &a1);
	    a2 = a * fhigh;
	    if (a2 > upl) {
		a2 = upl / a;
		if (a2 > flow)
		    vectd_(&a2, &upl);
	    } else {
		if (a2 > botl)
		    vectd_(&fhigh, &a2);
	    }
	}
	setusv_("LW", &two_thousand, 2);
	gcolor_(&ione);
	a = fhigh * 0.9;
	upl *= 0.5;
	cntr = 1.;
	size = 12.;
	angd = 0.;
	(void) sprintf(label, "VARIANCE=%15.7g", variance);
	printf(" FFT total variance=%12.4g, <20-km=%12.4g\n", variance,
	       variance_20km);
	label_length = 24;
	/* plchhq_(&a, &upl, label, &size, &angd, &cntr, label_length); */
	head_(&date, &start_time, &end_time, &dummy, 1);
	titleline = 1;
	fnote_(&titleline, var, strlen(var));
	idline_(&titleline, label, strlen(label));
	titleline = 2;
	fnote_(&titleline, "FFT", 3);
	sprintf(label, "pts/seg=%6d, window=%2d, smooth bins=%5d",
		segment_length, iwindow + 1, smooth_bins);
	idline_(&titleline, label, 44);
	if (!strncmp(var, "WIND", 4)) {
	    titleline = 3;
	    sprintf(label, "wind in %d degree direction", proj_d);
	    idline_(&titleline, label, strlen(label));
	}
	gcolor_(&ione);

	/*
	 * this strange call forces color reset, needed if plot frames are
	 * overlaid by 'med'
	 */
	line_(&thousand, &thousand, &thousand, &thousand);
	frame_();
    }
    if (ish & 8) {
	/* set limits */
	upl = 1.e-10;
	botl = 1.e10;
	for (isgn = -1; isgn <= 1; isgn++) {
	    if (isgn && !show_errors)
		continue;
	    if (sp_var == 0 || !strncmp (variable, "UX", 2))
		ae = 2.0;
	    else
		ae = 1.5;
//	    ae = pow(ae, (double) 1.5) * TWOPI / tas_average;
	    ae *= pow ((double)(TWOPI/tas_average), TWOTHIRDS);

	    /* now average together values within smooth_bins */
	    df = (log10((double) fhigh) - log10((double) flow)) / smooth_bins;
	    logf = log10((double) nu[2]);
//	    a = pow((double) ps[1] * (1. + isgn / sqrt(sigma)), (double) (1.5))
//		* ae * pow((double) nu[1], (double) 2.5);
//		correction Oct 2013 (but still need to fix the error limits)
	    a = ps[1] * ae * pow ((double) nu[1], (double) FIVETHIRDS);
	    a = pow (a, 1.5);
	    if (botl > a)
		botl = a;
	    if (upl < a)
		upl = a;
	    i = 1;
	    a = 0.;
	    frequency = 0.;
	    m = 0;
	    while (i++ < max2) {
		while (log10((double) nu[i]) - logf > df) {
		    if (m != 0) {
			frequency /= m;
//			a = pow((double) a * (1. + isgn / sqrt(sigma * m)) / m,
//				(double) (1.5))
//			    * ae * pow((double) frequency, (double) 2.5);
//				correction Oct 2013
			a /= m;
			a = pow (a, 1.5);
			if (botl > a)
			    botl = a;
			if (upl < a)
			    upl = a;
			frequency = a = 0.;
			m = 0;
		    }
		    logf += df;
		}
		frequency += nu[i];
		a += ps[i] * ae * pow ((double) nu[i], (double) FIVETHIRDS);
		m++;

		if (smooth_bins == MAXSMOOTH) {

		    /*
		     * convert to equiv of eddy diss. rate
		     */
//		    a = pow((double) a * (1. + isgn / sqrt(sigma)),
//			    (double) (1.5))
//			* ae * pow((double) frequency, (double) 2.5);
//			correction Oct 2013
		    a = pow (a, 1.5);
		    if (botl > a)
			botl = a;
		    if (upl < a)
			upl = a;
		    frequency = a = 0.;
		    m = 0;
		}
	    }
	    if (m != 0) {
		frequency /= m;
//		a = pow((double) a * (1. + isgn / sqrt(sigma * m)) / m,
//			(double) (1.5))
//		    * ae * pow((double) frequency, (double) 2.5);
			// changed Oct 2013
		a = pow ((double) (a/m), 1.5);
		if (botl > a)
		    botl = a;
		if (upl < a)
		    upl = a;
		frequency = a = 0.;
		m = 0;
	    }
	}
	if (botl < 1.e-12 && 1.e-12 < upl)
	    botl = 1.e-12;
	a = log10((double) botl);
	if (a < 0)
	    i = a - 1.;
	else
	    i = a;
	botl = pow((double) 10., (double) i);
	if (upl > 1.e6 && 1.e6 > botl)
	    upl = 1.e6;
	a = log10((double) upl);
	if (a > 0)
	    i = a + 1.;
	else
	    i = a;
	upl = pow((double) 10., (double) i);
	if (epslow_spec > 0.)
	    botl = epslow_spec;
	if (epshigh_spec > 0.)
	    upl = epshigh_spec;
	ncar_(&flow, &fhigh, &botl, &upl, &ione, &ione, &ione, &ione, &ifour, &ione);
	lablx_(":PGL:T:RU: :L:(:U:H:L:Z):U:", 27);
	lably_("E(:GL:T:RU:) [m:S1:2 s:S2:-3  ]", 31);
	add_lscale(&flow, &fhigh, &botl, &upl);
	for (isgn = -1; isgn <= 1; isgn++) {
	    if (isgn && !show_errors)
		continue;
	    if (isgn) {
		gcolor_(&ithree);
		setusv_("LW", &thousand, 2);
	    } else {
		gcolor_(&icolor);
		setusv_("LW", &two_thousand, 2);
	    }
	    if (sp_var == 0 || !strncmp (variable, "UX", 2))
		ae = 2.0;
	    else
		ae = 1.5;
//	    ae = pow(ae, (double) 1.5) * TWOPI / tas_average;
//	    		correction Oct 2013
	    ae *= pow ((double)(TWOPI/tas_average), (double) TWOTHIRDS);
	    /* now average together values within smooth_bins */
	    df = (log10((double) fhigh) - log10((double) flow)) / smooth_bins;
	    logf = log10((double) nu[1]);
//	    a = pow((double) ps[1] * (1. + isgn / sqrt(sigma)), (double) (1.5))
//		* ae * pow((double) nu[1], (double) 2.5);
//			correction Oct 2013
	    a = ps[1] * ae * pow ((double) nu[1], (double) FIVETHIRDS);
	    a = pow (a, 1.5);
	    if (a < botl)
		a = botl;
	    frstpt_(nu + 1, &a);
	    i = 1;
	    a = 0.;
	    frequency = 0.;
	    m = 0;
	    while (i++ < max2) {
		while (log10((double) nu[i]) - logf > df) {
		    if (m != 0) {
			frequency /= m;
//			a = pow((double) a * (1. + isgn / sqrt(sigma * m)) / m,
//				(double) (1.5))
//			    * ae * pow((double) frequency, (double) 2.5);
//			        correction Oct 2013
			a = pow ((double) (a/m), (double) 1.5);
			if (a < botl)
			    a = botl;
			vector_(&frequency, &a);
			frequency = a = 0.;
			m = 0;
		    }
		    logf += df;
		}
		frequency += nu[i];
	        a += ps[i] * ae * pow ((double) nu[i], (double) FIVETHIRDS);
		m++;
		if (smooth_bins == MAXSMOOTH) {

		    /*
		     * convert to equiv of eddy diss. rate
		     */
/*		    a = pow((double) a * (1. + isgn / sqrt(sigma)),
			    (double) (1.5))
			* ae * pow((double) frequency, (double) 2.5);
 */
		    a = pow (a, 1.5);
		    if (a < botl)
			a = botl;
		    vector_(&frequency, &a);
		    frequency = a = 0.;
		    m = 0;
		}
	    }
	    if (m != 0) {
		frequency /= m;
/*
		a = pow((double) a * (1. + isgn / sqrt(sigma * m)) / m,
			(double) (1.5))
		    * ae * pow((double) frequency, (double) 2.5);
 */
		a = pow (a/m, 1.5);
		if (a < botl)
		    a = botl;
		vector_(&frequency, &a);
		frequency = a = 0.;
		m = 0;
	    }
	}
	setusv_("LW", &two_thousand, 2);
	gcolor_(&iseven);
	frstpt_(nu + max2, &eps_bar);
	a = nu[max2] / 10.;
	vector_(&a, &eps_bar);
	gcolor_(&ione);
	head_(&date, &start_time, &end_time, &dummy, 1);
	titleline = 1;
	fnote_(&titleline, var, strlen(var));
	titleline = 2;
	fnote_(&titleline, "FFT", 3);
	sprintf(label, ":GL:E:RU:=%.3e :GL:,:RU:%.3e", eps_bar, eps2_bar);
	idline_(&titleline, label, strlen(label));
	if (!strncmp(var, "WIND", 4)) {
	    titleline = 3;
	    sprintf(label, "wind in %d degree direction", proj_d);
	    idline_(&titleline, label, strlen(label));
	}
	gcolor_(&ione);

	/*
	 * this strange call forces color reset, needed if plot frames are
	 * overlaid by 'med'
	 */
	line_(&thousand, &thousand, &thousand, &thousand);
	frame_();
    }
    free(ps);
    free(nu);
    return;
}

void
cosplt(dps1, dps2, sigma, sp_var1, var1, sp_var2, var2,
       proj_d1, proj_d2)
    double          dps1[], dps2[];
    double          sigma;
    int             sp_var1, sp_var2;
    char            var1[], var2[];
    int             proj_d1, proj_d2;

{

    float          *nu, *ps, frequency;
    float           upl, botl;
    float           ps_plot;
    float           flow = .001;
    float           fhigh = 10.;
    float           a, a1, a2;
    float           ae;
    int             max2, m, i, j;
    int             label_length;
    int             ipc, ipcold = 60;
    int             isgn;
    int             m_eps = 0;
    float           eps, eps_bar = 0., eps2_bar = 0.;
    char            label[50];
    float           df, logf;
    int             iset;

    max2 = segment_length / 2;
    nu = vector(1, max2);	/* allocate storage */
    ps = vector(1, max2);
    /* plot cospectrum and quadrature */
    for (iset = 0; iset < 2; iset++) {
	upl = -1.e10, botl = 1.e10;
	flux = 0.;
	for (i = 1; i <= max2; i++) {
	    nu[i] = i * ((float) sample_frequency) / segment_length;
	    if (nu[i] > flux_min_freq && nu[i] <= flux_max_freq)
	        flux += dps1[i];
// should include pad_weight here as for main spectrum -- do later
	    if (iset == 0)
		ps[i] = (float) nu[i] * (dps1[i] * (((double) segment_length) / sample_frequency));
	    else
		ps[i] = (float) nu[i] * (dps2[i] * (((double) segment_length) / sample_frequency));
	    if (upl < ps[i])
		upl = ps[i];
	    if (botl > ps[i])
		botl = ps[i];

	    ipc = i * 10 / max2 + 60;
	    if (ipc != ipcold) {
		ipcold = ipc;
	    }
	}
	fhigh = 10.;
	if (sample_frequency > 10.) {
	    fhigh *= 10.;
	}
	if (flow_spec > 0.)
	    flow = flow_spec;
	if (fhigh_spec > 0.)
	    fhigh = fhigh_spec;
	set_limits(&botl, &upl);
	if (cplow_spec != 0.)
	    botl = cplow_spec;
	if (cphigh_spec != 0.)
	    upl = cphigh_spec;
	ncar_(&flow, &fhigh, &botl, &upl, &ione, &ione, &ifour, &itwo, &ithree, &ione);
	lablx_(":PGL:T:RU: :L:(:U:H:L:Z):U:", 27);
/*	lably_("P(:GL:T:RU:) [m:S1:2 s:S2:-1  ]", 31); */
	ilbl = 0;
	    if (!strncmp(var1, "WI", 2)) ilbl++;
	    if (!strncmp(var1, "LONW", 4)) ilbl++;
	    if (!strncmp(var1, "LATW", 4)) ilbl++;
	    if (!strncmp(var2, "WI", 2)) ilbl++;
	    if (!strncmp(var2, "LONW", 4)) ilbl++;
	    if (!strncmp(var2, "LATW", 4)) ilbl++;
	    if (ilbl < 2) ilbl = 0;
	    if (ilbl) {
	        if (iset == 0)
	            lably_(":GL:T:RU:C(:GL:T:RU:) [m:S1:2 s:S2:-2  ]", 40);
	        else
	            lably_(":GL:T:RU:Q(:GL:T:RU:) [m:S1:2 s:S2:-2  ]", 40);
	    } else {
	        if (iset == 0)
	            lably_(":GL:T:RU:C(:GL:T:RU:)          ", 31);
	        else
	            lably_(":GL:T:RU:Q(:GL:T:RU:)          ", 31);
	    }
	add_lscale(&flow, &fhigh, &botl, &upl);
	df = (log10((double) fhigh) - log10((double) flow)) / smooth_bins;
	for (isgn = -1; isgn <= 1; isgn++) {
	    if (isgn && !show_errors)
		continue;
	    if (isgn) {
		gcolor_(&ithree);
		setusv_("LW", &thousand, 2);
	    } else {
		gcolor_(&icolor);
		setusv_("LW", &two_thousand, 2);
	    }
	    /* now average together values within smooth_bins */
	    logf = log10((double) nu[2]);
	    a = ps[1] * (1. + isgn / sqrt(sigma));
	    if (a < botl)
		a = botl;
	    frstpt_(nu + 1, &a);
	    i = 1;
	    a = 0.;
	    frequency = 0.;
	    m = 0;
	    while (i++ < max2) {
		while (log10((double) nu[i]) - logf > df) {
		    if (m) {
			a *= (1. + isgn / sqrt(sigma * m)) / m;
			frequency /= m;
			if (a < botl)
			    a = botl;
			vector_(&frequency, &a);
			a = frequency = 0.;
			m = 0;
		    }
		    logf += df;
		}
		frequency += nu[i];
		a += ps[i];
		m++;
		if (smooth_bins == MAXSMOOTH) {
		    a *= (1. + isgn / sqrt(sigma));
		    vector_(&frequency, &a);
		    frequency = a = 0.;
		    m = 0;
		}
	    }
	    if (m) {
		a *= (1. + isgn / sqrt(sigma * m)) / m;
		frequency /= m;
		if (a < botl)
		    a = botl;
		vector_(&frequency, &a);
	    }
	}
	setusv_("LW", &two_thousand, 2);
	gcolor_(&ione);
	a = fhigh * 0.9;
	printf("FLUX=%15.7g, limits %f - %f Hz\n", 
		flux, flux_min_freq, flux_max_freq);
	(void) sprintf(label, "FLUX=%15.7g, limits %f - %f Hz", 
		flux, flux_min_freq, flux_max_freq);
	titleline = 1;
	idline_(&titleline, label, strlen(label));
	head_(&date, &start_time, &end_time, &dummy, 1);
	sprintf(label, "%s x %s", var1, var2);
	fnote_(&titleline, label, strlen(label));
	titleline = 2;
	sprintf(label, "FFT, %s delayed %d ms", var1, delay);
	fnote_(&titleline, label, strlen(label));
	sprintf(label, "pts/seg=%6d, window=%2d, smooth bins=%5d",
		segment_length, iwindow + 1, smooth_bins);
	idline_(&titleline, label, 44);
	if (!strncmp(var1, "WIND", 4)) {
	    titleline = 3;
	    sprintf(label, "wind in %d degree direction", proj_d1);
	    idline_(&titleline, label, strlen(label));
	}
	if (!strncmp(var2, "WIND", 4)) {
	    titleline = 3;
	    sprintf(label, "wind in %d degree direction", proj_d2);
	    idline_(&titleline, label, strlen(label));
	}
	gcolor_(&ione);

	/*
	 * this strange call forces color reset, needed if plot frames are
	 * overlaid by 'med'
	 */
	line_(&thousand, &thousand, &thousand, &thousand);
	frame_();
    }
    free(ps);
    free(nu);
    return;
}

void
cophase(dps1, dps2, vps1, vps2, sigma, sp_var1, var1, sp_var2, var2,
	proj_d1, proj_d2)
    double          dps1[], dps2[];
    double          vps1[], vps2[];
    double          sigma;
    int             sp_var1, sp_var2;
    char            var1[], var2[];
    int             proj_d1, proj_d2;

{

    float          *nu, *ps, frequency;
    float           upl = -1.e10, botl = 1.e10;
    float           ps_plot;
    float           flow = .001;
    float           fhigh = 10.;
    float           a, a1, a2;
    double          c1, c2, c3, c4;
    double          b;
    float           ae;
    int             max2, m, i, j;
    int             label_length;
    int             ipc, ipcold = 60;
    int             isgn;
    int             m_eps = 0;
    float           eps, eps_bar = 0., eps2_bar = 0.;
    char            label[50];
    float           df, logf;
    int             iset, ifrst;

    printf(" entered cophase subroutine \n");
    max2 = segment_length / 2;
    nu = vector(1, max2);	/* allocate storage */
    ps = vector(1, max2);
    /* plot coherence and phase */
    /* (note: +ve phase means cvar is ahead of var) */
    for (iset = 0; iset < 2; iset++) {
	for (i = 1; i <= max2; i++) {
	    nu[i] = i * ((float) sample_frequency) / segment_length;
	}
	fhigh = 10.;
	if (sample_frequency > 10.) {
	    fhigh *= 10.;
	}
	if (flow_spec > 0.)
	    flow = flow_spec;
	if (fhigh_spec > 0.)
	    fhigh = fhigh_spec;
	if (iset == 0) {
	    botl = 0.;
	    upl = 1.;
	} else {
	    botl = phlow_spec;
	    upl = phhigh_spec;
	}
	ncar_(&flow, &fhigh, &botl, &upl, &ione, &ione, &ifour, &itwo, &ithree, &ione);
	lablx_(":PGL:T:RU: :L:(:U:H:L:Z):U:", 27);
/*	lably_("P(:GL:T:RU:) [m:S1:2 s:S2:-1  ]", 31); */
	if (iset == 0)
	    lably_("Coh(:GL:T:RU:)", 14);
	else
	    lably_(":G:U:R:(:GL:T:RU:)", 18);
	add_lscale(&flow, &fhigh, &botl, &upl);
	df = (log10((double) fhigh) - log10((double) flow)) / smooth_bins;
	for (isgn = -1; isgn <= 1; isgn++) {
	    if (isgn && !show_errors)
		continue;
	    if (isgn) {
		gcolor_(&ithree);
		setusv_("LW", &thousand, 2);
	    } else {
		gcolor_(&icolor);
		setusv_("LW", &two_thousand, 2);
	    }
	    /* now average together values within smooth_bins */
	    logf = log10((double) nu[1]);
	    ifrst = 1;
	    i = 1;
	    a = 0.;
	    c1 = c2 = c3 = c4 = 0.;
	    frequency = 0.;
	    m = 0;
	    while (i++ < max2) {
		while (log10((double) nu[i]) - logf > df) {
		    if (m) {
			if (iset == 0) {
			    a = (float) sqrt((c1 * c1 + c2 * c2) / (c3 * c4));
			    a += (isgn * sqrt((1. - a * a) / (2. * sigma * m)));
			} else {
			    b = isgn * 0.7071068 / sqrt(sigma * m);
			    a = (float) (180. / PI
				     * atan2(c2 * (1. + b), c1 / (1. + b)));
			    while (a > upl) a -= 360.;
			    while (a < botl) a += 360.;
			}
			frequency /= m;
			if (a < botl)
			    a = botl;
			if (ifrst) {
			    ifrst = 0;
			    frstpt_(&frequency, &a);
			    if (iset == 0) {
				printf ("%lf %lf\n", frequency, a);
			    }
			} else
			    vector_(&frequency, &a);
			    if (iset == 0) {
				printf ("%lf %lf\n", frequency, a);
			    }
			    if (iset == 1) {
//		    	        printf("%lf %lf\n", frequency, a);
			    }
			a = frequency = 0.;
			c1 = c2 = c3 = c4 = 0.;
			m = 0;
		    }
		    logf += df;
		}
		frequency += nu[i];
		c1 += dps1[i];
		c2 += dps2[i];
		c3 += vps1[i];
		c4 += vps2[i];
		m++;
		if (smooth_bins == MAXSMOOTH) {
		    if (iset == 0) {
			a = (float) sqrt((c1 * c1 + c2 * c2)
					 / (c3 * c4));
			a += (isgn * sqrt((1. - a * a) / (2. * sigma)));
		    } else {
			b = isgn * 0.7071068 / sqrt(sigma);
			a = (float) (180. / PI
				     * atan2(c2 * (1. + b), c1 / (1. + b)));
			    while (a > upl) a -= 360.;
			    while (a < botl) a += 360.;
		    }
		    if (ifrst) {
			ifrst = 0;
			frstpt_(&frequency, &a);
		    } else
			vector_(&frequency, &a);
			    if (iset == 1) {
		    	        printf("%lf %lf\n", frequency, a);
			    }
		    frequency = a = 0.;
		    c1 = c2 = c3 = c4 = 0.;
		    m = 0;
		}
	    }
	    if (m) {
		if (iset == 0) {
		    a = (float) sqrt((c1 * c1 + c2 * c2)
				     / (c3 * c4));
		    a += (isgn * sqrt((1. - a * a) / (2. * sigma * m)));
		} else {
		    b = isgn * 0.7071068 / sqrt(sigma * m);
		    a = (float) (180. / PI
				 * atan2(c2 * (1. + b), c1 / (1. + b)));
			    while (a > upl) a -= 360.;
			    while (a < botl) a += 360.;
		}
		frequency /= m;
		if (a < botl)
		    a = botl;
		vector_(&frequency, &a);
			    if (iset == 1) {
		    	        printf("%lf %lf\n", frequency, a);
			    }
	    }
	}
	setusv_("LW", &two_thousand, 2);
	gcolor_(&ione);
	a = fhigh * 0.9;
	titleline = 1;
	head_(&date, &start_time, &end_time, &dummy, 1);
	(void) sprintf(label, "%s delayed by %d ms", 
		var1, delay);
	titleline = 1;
	idline_(&titleline, label, strlen(label));
	sprintf(label, "%s x %s", var1, var2);
	fnote_(&titleline, label, strlen(label));
	titleline = 2;
	fnote_(&titleline, "FFT", 3);
	sprintf(label, "pts/seg=%6d, window=%2d, smooth bins=%5d",
		segment_length, iwindow + 1, smooth_bins);
	idline_(&titleline, label, 44);
	if (!strncmp(var1, "WIND", 4)) {
	    titleline = 3;
	    sprintf(label, "wind in %d degree direction", proj_d1);
	    idline_(&titleline, label, strlen(label));
	}
	if (!strncmp(var2, "WIND", 4)) {
	    titleline = 3;
	    sprintf(label, "wind in %d degree direction", proj_d2);
	    idline_(&titleline, label, strlen(label));
	}
	gcolor_(&ione);

	/*
	 * this strange call forces color reset, needed if plot frames are
	 * overlaid by 'med'
	 */
	line_(&thousand, &thousand, &thousand, &thousand);
	frame_();
    }
    free(ps);
    free(nu);
    return;
}

void
cphase(nu, dps1, dps2, vps1, vps2, sigma, sp_var1, var1, sp_var2, var2,
       proj_d1, proj_d2)
    float           nu[];
    float           dps1[], dps2[];
    float           vps1[], vps2[];
    double          sigma;
    int             sp_var1, sp_var2;
    char            var1[], var2[];
    int             proj_d1, proj_d2;

{

    float          *ps, frequency;
    float           upl = -1.e10, botl = 1.e10;
    float           ps_plot;
    float           flow = .001;
    float           fhigh = 10.;
    float           a, a1, a2;
    double          c1, c2, c3, c4;
    double          b;
    float           ae;
    int             max2, m, i, j;
    int             label_length;
    int             ipc, ipcold = 60;
    int             isgn, ifrst;
    int             m_eps = 0;
    float           eps, eps_bar = 0., eps2_bar = 0.;
    char            label[50];
    float           df, logf;
    int             iset;

    max2 = 1. / resolution;
    ps = vector(1, max2);
    /* plot coherence and phase */
    for (iset = 0; iset < 2; iset++) {
	flow = 0.001, fhigh = 10.;
	if (sample_frequency > 10.) {
	    fhigh *= 10.;
	}
	if (flow_spec > 0.)
	    flow = flow_spec;
	if (fhigh_spec > 0.)
	    fhigh = fhigh_spec;
	if (iset == 0) {
	    botl = 0.;
	    upl = 1.;
	} else {
	    botl = -180.;
	    upl = 180.;
	}
	ncar_(&flow, &fhigh, &botl, &upl, &ione, &ione, &ifour, &itwo, &ithree, &ione);
	lablx_(":PGL:T:RU: :L:(:U:H:L:Z):U:", 27);
/*	lably_("P(:GL:T:RU:) [m:S1:2 s:S2:-1  ]", 31); */
	if (iset == 0)
	    lably_("Coh(:GL:T:RU:)", 14);
	else
	    lably_(":G:U:R:(:GL:T:RU:)", 18);
	add_lscale(&flow, &fhigh, &botl, &upl);
	df = (log10((double) fhigh) - log10((double) flow)) / smooth_bins;
	for (isgn = -1; isgn <= 1; isgn++) {
	    if (isgn) continue;			/* suppress error limits; FIX */
	    if (isgn && !show_errors)
		continue;
	    if (isgn) {
		gcolor_(&ithree);
		setusv_("LW", &thousand, 2);
	    } else {
		gcolor_(&icolor);
		setusv_("LW", &two_thousand, 2);
	    }
	    /* now average together values within smooth_bins */
	    logf = log10((double) nu[2]);
	    ifrst = 1;
	    i = 1;
	    a = 0.;
	    c1 = c2 = c3 = c4 = 0.;
	    frequency = 0.;
	    m = 0;
	    while (i++ < max2) {
		while (log10((double) nu[i]) - logf > df) {
		    if (m) {
			if (iset == 0) {
			    a = (float) sqrt((c1 * c1 + c2 * c2) / (c3 * c4));
			    a += (isgn / sigma * sqrt((double) (1. - a * a) / (2. * m)));
			} else {
			    b = isgn * 0.7071068 / (sigma * sqrt((double) m));
			    a = (float) (180. / PI
				     * atan2(c2 * (1. + b), c1 / (1. + b)));
			}
			frequency /= m;
			if (a < botl)
			    a = botl;
			if (ifrst) {
			    frstpt_(&frequency, &a);
			    ifrst = 0;
			} else
			    vector_(&frequency, &a);
			a = frequency = 0.;
			c1 = c2 = c3 = c4 = 0.;
			m = 0;
		    }
		    logf += df;
		}
		frequency += nu[i];
		c1 += dps1[i];
		c2 += dps2[i];
		c3 += vps1[i];
		c4 += vps2[i];
		m++;
		if (smooth_bins == MAXSMOOTH) {
		    if (iset == 0) {
			a = (float) sqrt((c1 * c1 + c2 * c2) / (c3 * c4));
			if (a < 1.) {
			    a += (isgn / sigma) * sqrt((1. - a * a) / 2.);
			}
		    } else {
			b = isgn * 0.7071068 / sigma;
			a = (float) (180. / PI
				     * atan2(c2 * (1. + b), c1 / (1. + b)));
		    }
		    if (ifrst) {
			frstpt_(&frequency, &a);
			ifrst = 0;
		    } else
			vector_(&frequency, &a);
		    frequency = a = 0.;
		    c1 = c2 = c3 = c4 = 0.;
		    m = 0;
		}
	    }
	    if (m) {
		if (iset == 0) {
		    a = (float) sqrt((c1 * c1 + c2 * c2)
				     / (c3 * c4));
		    a += (isgn / sigma) * sqrt((1. - a * a) / (2. * m));
		} else {
		    b = isgn * 0.7071068 / ((sigma * sqrt((double) m)) * sqrt((double) m));
		    a = (float) (180. / PI
				 * atan2(c2 * (1. + b), c1 / (1. + b)));
		}
		frequency /= m;
		if (a < botl)
		    a = botl;
		vector_(&frequency, &a);
	    }
	}
	setusv_("LW", &two_thousand, 2);
	gcolor_(&ione);
	a = fhigh * 0.9;
	titleline = 1;
	head_(&date, &start_time, &end_time, &dummy, 1);
	(void) sprintf(label, "%s delayed by %d ms", 
		var1, delay);
	titleline = 1;
	idline_(&titleline, label, strlen(label));
	sprintf(label, "%s x %s", var1, var2);
	fnote_(&titleline, label, strlen(label));
	titleline = 2;
	fnote_(&titleline, "MEM", 3);
	(void) sprintf(label, "poles=%5d, resln=%8.6f, smooth bins=%5d",
		       poles, resolution, smooth_bins);
	idline_(&titleline, label, 46);
	if (!strncmp(var1, "WIND", 4)) {
	    titleline = 3;
	    sprintf(label, "wind in %d degree direction", proj_d1);
	    idline_(&titleline, label, strlen(label));
	}
	if (!strncmp(var2, "WIND", 4)) {
	    titleline = 3;
	    sprintf(label, "wind in %d degree direction", proj_d2);
	    idline_(&titleline, label, strlen(label));
	}
	gcolor_(&ione);

	/*
	 * this strange call forces color reset, needed if plot frames are
	 * overlaid by 'med'
	 */
	line_(&thousand, &thousand, &thousand, &thousand);
	frame_();
    }
    free(ps);
    return;
}

void
set_limits(botl, upl)
    float          *botl, *upl;

{
    float           diff, dtest, low, high;
    int             id;

    printf(" entry to set_limits with limits %f to %f\n", *botl, *upl);
    diff = *upl - *botl;
    dtest = log10((double) diff);
    if (dtest > 0)
	id = dtest + 1.;
    else
	id = dtest;
    dtest = pow((double) 10., (double) id);
    low = *botl / dtest;
    high = *upl / dtest;
    if (low < -0.5)
	low = -1.;
    else if (low < -0.2)
	low = -0.5;
    else
	low = -0.2;
    if (high > 0.5)
	high = 1.;
    else if (high > 0.2)
	high = 0.5;
    else
	high = 0.2;
    if (high > -low)
	low = -high;
    if (low < -high)
	high = -low;
    *botl = low * dtest;
    *upl = high * dtest;
    printf(" limits set to %f and %f, dtest = %f, id = %d\n", *botl, *upl,
	   dtest, id);
    return;
}

float
ftsec(t)
    float           t;

{
    int             itim;
    float           ans;

    itim = t;
    ans = (itim % 100 + ((itim / 100) % 100) * 60 + (itim / 10000) * 3600) + (t - itim);
    return (ans);
}

float
fsect(t)
    float	    t;
{
    float	    f;
    int	  	    it;

    it = t;
    f = t - it;
    return (f + it % 60 + ((it / 60) % 60) * 100 + (it / 3600) * 10000);
}

int
itsec(t)
    int             t;

{
    return (t % 100 + ((t / 100) % 100) * 60 + (t / 10000) * 3600);
}

int
isect(ts)
    int             ts;

{
    return (ts % 60 + ((ts / 60) % 60) * 100 + (ts / 3600) * 10000);
}


void
detrend(data_b, mean_local, trend_local, set_zero)
    float          *data_b;
    double          mean_local, trend_local;
    int		    set_zero;

{
    int             i, nRMS;
    float	    midpoint;

//  midpoint = ((duration - 1.) / (0.002 * dhd1.ideltt));
    midpoint = ((duration) / (0.002 * dhd1.ideltt));
    LogRMS = 0.;
    nRMS = 0;
    for (i = 0; i < points; i++) {
	if( *(data_b + i) != MISSING_DATA) {
	    *(data_b + i) -= (mean_local + trend_local * (i - midpoint));
	    LogRMS += *(data_b+i) * *(data_b+i);
	    nRMS++;
	} else if (set_zero) { *(data_b + i) = 0.;}
    }
    LogRMS /= nRMS+1.e-5;
    LogRMS = pow(LogRMS, 0.5);
    return;
}

void
add_lscale(flow, fhigh, botl, upl)
    float          *flow, *fhigh, *botl, *upl;

{
    int             i, j;
    double          a, b;
    float           f, y1, y2, y3, y4;
    char            label[30];
    float           fl, fr, fb, ft, ul, ur, ub, ut;
    int             ll;
    float           angd = 0., cntr = 0., size = 10.;	/* arguments for plchhq */

    y1 = *upl;
    getset_(&fl, &fr, &fb, &ft, &ul, &ur, &ub, &ut, &ll);
    if (ll == 4) {
	a = (log10((double) *upl) - log10((double) *botl)) / 60.;
	y2 = y1 * pow(10., a);
	a *= 1.4;
	y3 = y1 * pow(10., a);
	a *= 1.5;
	y4 = y1 * pow(10., a);
    } else {
	a = (*upl - *botl) / 60.;
	y2 = y1 + a;
	y3 = y1 + 1.4 * a;
	y4 = y1 + (1.4 * 1.5) * a;
    }
    for (i = -1; i <= 1; i++) {
	b = i;
	a = pow(10., b);
	f = tas_average / (a * 1000.);
	line_(&f, &y1, &f, &y3);
	if (i < 1)
	    sprintf(label, "%3.1f", a);
	else
	    sprintf(label, "%2.0f", a);
	plchhq_(&f, &y4, label, &size, &angd, &cntr, strlen(label));
	if (i == 1)
	    continue;
	for (j = 2; j < 10; j++) {
	    b = a * j;
	    f = tas_average / (b * 1000.);
	    line_(&f, &y1, &f, &y2);
	}
    }
    sprintf(label, "wavelength [km]:PRL:0 ");
    f = *flow;
    cntr = 0.;
    size = 8.;
    plchhq_(&f, &y4, label, &size, &angd, &cntr, strlen(label));
    return;
}
