executables = xoptimize_change_points_gfort.exe xxsplit_data_gfort.exe xsegmentation_gfort.exe xneighbors_gfort.exe
FC     = gfortran
FFLAGS = -O0 -Wall -Werror=unused-parameter -Werror=unused-variable -Werror=unused-function -Wno-maybe-uninitialized -Wno-surprising -fbounds-check -static -g -fmodule-private
obj    = kind.o constants.o random.o change_point_util.o split_data.o xoptimize_change_points.o xxsplit_data.o segmentation.o xsegmentation.o xneighbors.o

all: $(executables)

# Compile .f90 to .o
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

xoptimize_change_points_gfort.exe: kind.o constants.o random.o change_point_util.o split_data.o xoptimize_change_points.o
	$(FC) -o xoptimize_change_points_gfort.exe kind.o constants.o random.o change_point_util.o split_data.o xoptimize_change_points.o $(FFLAGS)

xxsplit_data_gfort.exe: kind.o random.o change_point_util.o split_data.o xxsplit_data.o
	$(FC) -o xxsplit_data_gfort.exe kind.o random.o change_point_util.o split_data.o xxsplit_data.o $(FFLAGS)

xsegmentation_gfort.exe: kind.o constants.o random.o segmentation.o xsegmentation.o
	$(FC) -o xsegmentation_gfort.exe kind.o constants.o random.o segmentation.o xsegmentation.o $(FFLAGS)

xneighbors_gfort.exe: change_point_util.o xneighbors.o
	$(FC) -o xneighbors_gfort.exe change_point_util.o xneighbors.o $(FFLAGS)

run: $(executables)
	./xoptimize_change_points_gfort.exe
	./xxsplit_data_gfort.exe
	./xsegmentation_gfort.exe
	./xneighbors_gfort.exe

clean:
	rm -f $(executables) $(obj)

