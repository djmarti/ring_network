#include "conf_file.h"
/*

   Some parts were taken from inih, released under the New BSD license (see
   LICENSE.txt). Go to the project home page for more info:

   http://code.google.com/p/inih/ 

 */

static struct option long_opts_1stround[] = {
    {"config-file", required_argument, NULL, 'c'},
    {0, 0, 0, 0}
};

static struct option long_opts[] = {
    /* These options don't set a flag.
       We distinguish them by their indices. */
    {"config-file", required_argument, NULL, 'c'},
    {"inputs-file", required_argument, NULL, 'x'},
    {"output-rate-file", required_argument, NULL, 'o'},
    {"output-adaptation-file", required_argument, NULL, 'A'},
    {"output-facilitation-file", required_argument, NULL, 'F'},
    {"duration-simulation", required_argument, NULL, 'T'},
    {"time-step", required_argument, NULL, 'D'},
    {"excitatory-weight ", required_argument, NULL, 'e'},
    {"inhibitory-weight", required_argument, NULL, 'i'},
    {"excitatory-concentration", required_argument, NULL, 'f'},
    {"inhibitory-concentration", required_argument, NULL, 'j'},
    {"normalize", no_argument, NULL, 'n'},
    {"verbose", no_argument, NULL, 'v'},
    {"help", no_argument, NULL, 'h'},
    {0, 0, 0, 0}
};

int process_options(void *usercfg, int argc, char *argv[])
{
    struct Parameters *pconfig = (struct Parameters *) usercfg;

    opterr = 1;                 /* report errors (default) */
    optind = 1;                 /* start again */

    int c;

    int option_index = 0;       /* getopt_long stores the option index here. */

    while ((c = getopt_long(argc, argv, "nhvc:x:o:A:F:T:D:e:i:f:j:t:d:s:z:",
                            long_opts, &option_index)) != -1) {
        switch (c) {
        case 0:
            if (strcmp(long_opts[option_index].name, "help") == 0) {
                usage(0, argv[0]);
                return 0;
            } else if (strcmp(long_opts[option_index].name, "normalize") == 0) {
                pconfig->normalize = 1;
            } else if (strcmp(long_opts[option_index].name, "verbose") == 0) {
                pconfig->verbose = 1;
            } else if (strcmp(long_opts[option_index].name, "output-rate-file") == 0) {
                strcpy(pconfig->output_rate_file, optarg);
            } else if (strcmp(long_opts[option_index].name, "output-adaptation-file") == 0) {
                strcpy(pconfig->output_adaptation_file, optarg);
            } else if (strcmp(long_opts[option_index].name, "output-facilitation-file") == 0) {
                strcpy(pconfig->output_facilitation_file, optarg);
            } else if (strcmp(long_opts[option_index].name, "inputs-file") == 0) {
                strcpy(pconfig->extinputs_file, optarg);
            } else if (strcmp(long_opts[option_index].name,
                              "duration-simulation") == 0) {
                pconfig->T = atof(optarg);
            } else if (strcmp(long_opts[option_index].name, "time-step") == 0) {
                pconfig->dt = atof(optarg);
            } else if (strcmp
                       (long_opts[option_index].name,
                        "excitatory-weight") == 0) {
                pconfig->kernel.J_E = atof(optarg);
            } else if (strcmp
                       (long_opts[option_index].name,
                        "inhibitory-weight") == 0) {
                pconfig->kernel.J_I = atof(optarg);
            } else if (strcmp
                       (long_opts[option_index].name,
                        "excitatory-concentration") == 0) {
                pconfig->kernel.m_E = atof(optarg);
            } else if (strcmp
                       (long_opts[option_index].name,
                        "inhibitory-concentration") == 0) {
                pconfig->kernel.m_I = atof(optarg);
            }
            break;
        case 'h':
            usage(0, argv[0]);
            return 0;
        case 'c':
            break;              /* we took care of that before */
        case 'v':
            pconfig->verbose = 1;
            break;
        case 'n':
            pconfig->normalize = 1;
            break;
        case 'x':
            strcpy(pconfig->extinputs_file, optarg);
            break;
        case 'o':
            strcpy(pconfig->output_rate_file, optarg);
            break;
        case 'A':
            strcpy(pconfig->output_adaptation_file, optarg);
            break;
        case 'F':
            strcpy(pconfig->output_facilitation_file, optarg);
            break;
        case 'T':
            pconfig->T = atof(optarg);
            break;
        case 'D':
            pconfig->dt = atof(optarg);
            break;
        case 'e':
            pconfig->kernel.J_E = atof(optarg);
            break;
        case 'i':
            pconfig->kernel.J_I = atof(optarg);
            break;
        case 'f':
            pconfig->kernel.m_E = atof(optarg);
            break;
        case 'j':
            pconfig->kernel.m_I = atof(optarg);
            break;
        default:
            printf("Invalid option -- '%s'\n", optarg);
            printf("Try `ring_rate --help for more information\n");
            abort();
        }
    }
    if (argc - option_index != 2) {
        return 1;
    }
    return 0;
}

int process_config_in_options(void *usercfg, int argc, char *argv[])
/* abridged version of process_options needed to check config */
{
    struct Parameters *pconfig = (struct Parameters *) usercfg;

    int c;
    /* getopt_long stores the option index here. */
    int option_index = 0;
    /* Don't report errors. This is the first round */
    opterr = 0;

    while ((c = getopt_long(argc, argv, "c:",
                            long_opts_1stround, &option_index)) != -1) {
        switch (c) {
        case 0:
            if (strcmp(long_opts[option_index].name, "config-file") == 0) {
                if (conf_parse(optarg, handler, pconfig) < 0) {
                    printf("Can't load '%s'\n", optarg);
                    return 1;
                }
                strcpy(pconfig->conf_file, optarg);
            }
            break;
        case 'c':
            strcpy(pconfig->conf_file, optarg);
            break;
        }
    }
    return 0;
}

void usage(int status, char *s)
{
    if (status != 0) {
        fprintf(stderr, "Usage: %s [OPTION]...\n", s);
        fprintf(stderr, "Try `%s --help' for more information.\n", s);
    } else {
        printf("Usage: %s [OPTION]...\n", s);
        printf("\
\n\
Simulate continuous ring model.\n");
        printf("\
  -c, --config-file=FILE         read configuration parameters from FILE\n\
  -x, --inputs-file=FILE         read external input specification from FILE\n\
  -o, --output-file=FILE         save results in FILE\n");

        printf("\
\n\
Simulation:\n\
  -T, --duration-simulation=NUM  run simulation NUM time units\n\
  -D, --time-step=NUM            use Euler method with timestep NUM\n\
  -n, --normalize                normalize response to dynamic range in trial.\n");

        printf("\
\n\
Lateral connectivity:\n\
  -e, --excitatory-weight=NUM    set weight of the excitatory footprint\n\
  -i, --inhibitory-weight=NUM    set weight of the inhibitory footprint\n\
  -f, --excitatory-concentration=NUM  set concentration for the exc. footprint\n\
  -j, --inhibitory-concentration=NUM  set concentration for the inh. footprint\n");

        printf("\
\n\
Miscellaneous:\n\
  -v, --verbose                  be more verbose (show parameters)\n\
  -h, --help                     display this help and exit\n\n");
        exit(status);
    }
}
