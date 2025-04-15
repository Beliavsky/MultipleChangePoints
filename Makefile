executables = xneighbors_gfort.exe xoptimize_change_points_gfort.exe xsegmentation_gfort.exe xxsplit_data_gfort.exe
FC     = gfortran
FFLAGS = -O0 -Wall -Werror=unused-parameter -Werror=unused-variable -Werror=unused-function -Wno-maybe-uninitialized -Wno-surprising -fbounds-check -static -g -fmodule-private
obj    = change_point_util.o xneighbors.o kind.o constants.o random.o split_data.o xoptimize_change_points.o info_crit.o segmentation.o xsegmentation.o xxsplit_data.o

all: $(executables)

# Compile .f90 to .o
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

xneighbors_gfort.exe: change_point_util.o xneighbors.o
	$(FC) -o xneighbors_gfort.exe change_point_util.o xneighbors.o $(FFLAGS)

xoptimize_change_points_gfort.exe: change_point_util.o kind.o constants.o random.o split_data.o xoptimize_change_points.o
	$(FC) -o xoptimize_change_points_gfort.exe change_point_util.o kind.o constants.o random.o split_data.o xoptimize_change_points.o $(FFLAGS)

xsegmentation_gfort.exe: kind.o constants.o random.o info_crit.o segmentation.o xsegmentation.o
	$(FC) -o xsegmentation_gfort.exe kind.o constants.o random.o info_crit.o segmentation.o xsegmentation.o $(FFLAGS)

xxsplit_data_gfort.exe: change_point_util.o kind.o random.o split_data.o xxsplit_data.o
	$(FC) -o xxsplit_data_gfort.exe change_point_util.o kind.o random.o split_data.o xxsplit_data.o $(FFLAGS)

run: $(executables)
	./xneighbors_gfort.exe
	./xoptimize_change_points_gfort.exe
	./xsegmentation_gfort.exe
	./xxsplit_data_gfort.exe

clean:
	rm -f $(executables) $(obj)

