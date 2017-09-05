
void dump_memory_table();
void mymalloc_init();
void *mymalloc_fullinfo(const char *varname, size_t n, const char *func, const char *file, int line);
void *mymalloc_movable_fullinfo(void *ptr, const char *varname, size_t n, const char *func, const char *file, int line);
void myfree_fullinfo(void *p, const char *func, const char *file, int line);
void myfree_movable_fullinfo(void *p, const char *func, const char *file, int line);
void *myrealloc_fullinfo(void *p, size_t n, const char *func, const char *file, int line);
void *myrealloc_movable_fullinfo(void *p, size_t n, const char *func, const char *file, int line);

void init(int argc, char **argv);
void finish();

void read_snap();
void clear_snap();

void read_params();

void mpi_distribute_items_to_tasks(void *data, int task_offset, int *n_items, int *max_n, int item_size, int commtag);

void make_image();

void pspace();

void radial();

void halo();

void remove_part();

void cut_region();

void collapse();

int plot_block_present(int block);
void plot_block_label(int block, char *label);
double *plot_block_val(int block);

void box_center();
void box_rotate();
double get_H2_cool_rate();
double get_H2_esc_frac();
double get_3b_form_heat();
void compute_H2_opacity(double temp, double N_H2_eff, double *opac);

void tgchem_spline_eval(int nval, double *positions, double *values, int nnew, double *new_positions, double *new_values);
void tgchem_spline_coefficients(int nval, double *values, double *coefficients);
void tgchem_spline_derivatives(int nval, double *values, double *derivatives);
void endrun();
void mpi_printf(const char *fmt, ...);
int myflush(FILE *fstream);
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream);
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE * stream);
double mysort(void *base, size_t nel, size_t width, int (*compar) (const void *, const void *));
void terminate_args();
void terminate_block();
int imin(int a, int b);
int imax(int a, int b);
double dabs(double a);
double dmin(double a, double b);
double dmax(double a, double b);
double second();
double measure_time();
double timediff(double t0, double t1);
