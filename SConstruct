import os

srcs = Split("""ring_rate.c
                conf_file.c
                parser.c
                parameters.c
                ext_inputs.c
                colormap.c
                circ_buffer.c
                utils.c
                cubature-20101018/cubature.c""")

cflags = '-std=gnu99 -O3 --pedantic -W -Winline \
-Wstrict-prototypes -Wno-sign-conversion -Wshadow -Wpointer-arith -Wcast-qual \
-Wcast-align -Wwrite-strings -Wnested-externs \
-fshort-enums -fno-common'

opt = Environment(CFLAGS = cflags + ' -DHAVE_INLINE=1') #LINKFLAGS = ['-pg'])

# Libs
lpath = [os.path.join(os.environ['HOME'], 'lib'), '.']
opt.Library('ring_network', srcs)

libs = ['ring_network', 'm', 'gsl', 'gslcblas']
opt.Program('simulate_trajectory.c', LIBS=libs, LIBPATH=lpath)
#opt.Program('generate_colormap.c', LIBS=libs, LIBPATH=lpath)#
#opt.Program('simulate_many_trials.c', LIBS=libs, LIBPATH=lpath)
