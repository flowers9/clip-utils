# compiled under Solaris 8, using gcc 4.0.1, gmake 3.80, gnu-ld 2.16

OS := $(shell uname -s)

# use -m64 if the programs will need access to more than 4 GB of memory;
# otherwise, not using it will save memory

CXX     ?= g++
DEBUG    = -g -O2 #-DCOMPRESS_READS
# note that the uncompressed version of read.cc has diverged from the
# compressed version (has a few extra bits in it)
CPPFLAGS = $(extra_includes) -I./include $(DEBUG) -Wall -Wextra -Wpointer-arith -Wshadow \
	-Wundef -Wcast-qual -Woverloaded-virtual -Wsign-promo -Wuninitialized \
	-Wcast-align
	#-Wconversion

LDFLAGS  = $(extra_libs) $(DEBUG)

ifeq ($(OS), SunOS)
DEBUG += -mcpu=v9 -m64
# for compiling on hero (pants uses /bin/g++)
CXX = /opt/csw/bin/g++
#CPPFLAGS += -I/apps/include
#LDFLAGS += -R/apps/lib -R/apps/lib/sparcv9
endif

ifeq ($(OS), Linux)
DEBUG += -O3 -pthread -mfpmath=sse -march=native -flto -fno-fat-lto-objects -fno-builtin -mmmx -msse -msse2 -mssse3 -msse4.1 -msse4.2 -mpopcnt -mfxsr
ifeq ($(HOST), pc23-gsc)
# more recent pacbio libraries need newer c++, as well as includes and libraries
CPPFLAGS += -std=c++17
PACBIO_INCLUDES = -I./pbbam-2.1.0/include -I./pbcopper-2.0.0/include
# need to put -lpbcopper before -lpbbam, or it won't find the library at run time
PACBIO_LIBS = -L/home/smrtlink_beta/current/bundles/smrttools/current/private/pacbio/pbbam/lib -L/home/smrtlink_beta/current/bundles/smrttools/current/private/pacbio/pbcopper/lib -Wl,-R/home/smrtlink_beta/current/bundles/smrttools/current/private/pacbio/pbbam/lib -Wl,-R/home/smrtlink_beta/current/bundles/smrttools/current/private/pacbio/pbcopper/lib -lpbcopper -lpbbam
else
CPPFLAGS += -std=c++11
PACBIO_INCLUDES = -I/home/raid2/LINUXOPT/miniconda2a/include
PACBIO_LIBS = -L/home/raid2/LINUXOPT/miniconda2a/lib -Wl,-R/home/raid2/LINUXOPT/miniconda2a/lib -lpbbam -lhts
endif
endif

obj/%.o: %.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) -o $@ $<

.PHONY: all

all: bin/clip bin/histogram_hash bin/library_stats bin/mask_repeats_hash bin/qc_stats1 bin/qc_stats2 bin/targets bin/read_stats bin/read_histogram bin/phred_hist bin/parse_output bin/repair_sequence2 bin/compress_blat bin/repair_sequence3 bin/mask_repeats_hashn bin/histogram_hashn bin/check_barcodes bin/screen_blat bin/filter_blat bin/parse_output2 bin/screen_pairs bin/arachne_create_xml bin/extract_seq_and_qual bin/split_fasta bin/copy_dbs bin/print_hash bin/print_hashn bin/screen_reads bin/pacbio_read_stats bin/sort_blast bin/add_passes bin/find_kmers bin/add_quality bin/interleave bin/tee bin/chris_prep bin/kmer_matching_setup bin/kmer_matching bin/extract_bam_well_sizes bin/barcode_separation bin/filter_bam_alignments bin/split_bam bin/filter_bam bin/extract_good_read_names bin/dot_hash bin/dot_hashn bin/histogram_hashl bin/screen_kmers_by_ref bin/find_kmers_hashl bin/print_hashl bin/screen_kmers_by_lib

bin/chris_prep: obj/chris_prep.o obj/breakup_line.o obj/open_compressed.o obj/strtostr.o obj/write_fork.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

bin/tee: obj/tee.o obj/breakup_line.o obj/strtostr.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

bin/interleave: obj/interleave.o obj/open_compressed.o obj/breakup_line.o obj/strtostr.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

bin/screen_reads: obj/screen_reads.o obj/open_compressed.o obj/breakup_line.o obj/strtostr.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

bin/print_hash: obj/print_hash.o obj/hash.o obj/next_prime.o obj/open_compressed.o obj/write_fork.o obj/breakup_line.o obj/strtostr.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

bin/print_hashn: obj/print_hashn.o obj/get_name.o obj/hashn.o obj/hist_lib_hashn.o obj/next_prime.o obj/open_compressed.o obj/pattern.o obj/read.o obj/time_used.o obj/write_fork.o obj/breakup_line.o obj/strtostr.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

bin/copy_dbs: obj/copy_dbs.o obj/open_compressed.o obj/write_fork.o obj/breakup_line.o obj/strtostr.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

bin/screen_pairs: obj/screen_pairs.o obj/get_name.o obj/hashn.o obj/hist_lib_hashn.o obj/next_prime.o obj/open_compressed.o obj/pattern.o obj/read.o obj/read_lib.o obj/time_used.o obj/write_fork.o obj/breakup_line.o obj/strtostr.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

bin/filter_blat: obj/filter_blat.o obj/open_compressed.o obj/strtostr.o obj/breakup_line.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

bin/parse_output: obj/parse_output.o obj/open_compressed.o obj/strtostr.o obj/breakup_line.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

bin/parse_output2: obj/parse_output2.o obj/open_compressed.o obj/strtostr.o obj/breakup_line.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

bin/clip: obj/clip.o obj/breakup_line.o obj/get_name.o obj/open_compressed.o obj/pattern.o obj/read.o obj/read_file.o obj/strtostr.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

bin/histogram_hash: obj/open_compressed.o obj/get_name.o obj/hash.o obj/hist_lib_hash.o obj/histogram_hash.o obj/next_prime.o obj/pattern.o obj/read.o obj/read_file.o obj/time_used.o obj/write_fork.o obj/breakup_line.o obj/strtostr.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

bin/kmer_matching_setup: obj/kmer_matching_setup.o obj/get_name.o obj/hash.o obj/hash_read_hits.o obj/hist_lib_hash.o obj/kmer_lookup_info.o obj/next_prime.o obj/open_compressed.o obj/pattern.o obj/read.o obj/read_file.o obj/time_used.o obj/write_fork.o obj/breakup_line.o obj/strtostr.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

# specify version of readline because different installs have varying versions of the
# "current" version (from 6 to 8), but all have version 5
bin/kmer_matching: obj/kmer_matching.o obj/breakup_line.o obj/get_name.o obj/hash_read_hits.o obj/hist_lib_hash.o obj/kmer_lookup_info.o obj/next_prime.o obj/open_compressed.o obj/pattern.o obj/read.o obj/read_file.o obj/strtostr.o obj/time_used.o obj/write_fork.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS) -lreadline

bin/library_stats: obj/find_library.o obj/open_compressed.o obj/get_name.o obj/library_match.o obj/library_read_lib.o obj/library_stats.o obj/parse_read.o obj/pattern.o obj/pretty_print.o obj/read.o obj/read_lib.o obj/read_match.o obj/breakup_line.o obj/strtostr.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

bin/mask_repeats_hash: obj/breakup_line.o obj/open_compressed.o obj/get_name.o obj/hash.o obj/hist_lib_hash.o obj/mask_repeats_hash.o obj/next_prime.o obj/pattern.o obj/read.o obj/read_file.o obj/time_used.o obj/write_fork.o obj/breakup_line.o obj/strtostr.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

bin/dot_hash: obj/dot_hash.o obj/open_compressed.o obj/hash.o obj/next_prime.o obj/write_fork.o obj/breakup_line.o obj/strtostr.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

bin/dot_hashn: obj/dot_hashn.o obj/open_compressed.o obj/hashn.o obj/next_prime.o obj/write_fork.o obj/breakup_line.o obj/strtostr.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

bin/screen_kmers_by_ref: obj/screen_kmers_by_ref.o obj/open_compressed.o obj/hashl.o obj/hashl_metadata.o obj/next_prime.o obj/write_fork.o obj/breakup_line.o obj/strtostr.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

bin/screen_kmers_by_lib: obj/screen_kmers_by_lib.o obj/open_compressed.o obj/hashl.o obj/hashl_metadata.o obj/next_prime.o obj/write_fork.o obj/breakup_line.o obj/strtostr.o obj/time_used.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

bin/find_kmers_hashl: obj/find_kmers_hashl.o obj/open_compressed.o obj/hashl.o obj/hashl_metadata.o obj/next_prime.o obj/write_fork.o obj/breakup_line.o obj/strtostr.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

bin/print_hashl: obj/print_hashl.o obj/hashl.o obj/hashl_metadata.o obj/next_prime.o obj/open_compressed.o obj/breakup_line.o obj/strtostr.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

bin/qc_stats1: obj/open_compressed.o obj/get_name.o obj/pattern.o obj/pretty_print.o obj/qc_read.o obj/qc_read_lib.o obj/qc_stats1.o obj/read.o obj/read_lib.o obj/breakup_line.o obj/strtostr.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

bin/qc_stats2: obj/open_compressed.o obj/get_name.o obj/pattern.o obj/pretty_print.o obj/qc_stats2.o obj/read.o obj/read_lib.o obj/breakup_line.o obj/strtostr.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

bin/read_stats : obj/open_compressed.o obj/get_name.o obj/hash.o obj/hist_lib_hash.o obj/next_prime.o obj/pattern.o obj/read.o obj/read_lib.o obj/read_stats.o obj/time_used.o obj/write_fork.o obj/breakup_line.o obj/strtostr.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

bin/read_histogram : obj/open_compressed.o obj/get_name.o obj/hash.o obj/hist_lib_hash.o obj/next_prime.o obj/pattern.o obj/read.o obj/read_lib.o obj/read_histogram.o obj/time_used.o obj/write_fork.o obj/breakup_line.o obj/strtostr.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

bin/targets: obj/targets.o obj/open_compressed.o obj/get_name.o obj/pattern.o obj/read.o obj/read_lib.o obj/breakup_line.o obj/strtostr.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

bin/phred_hist: obj/open_compressed.o obj/phred_hist.o obj/pretty_print.o obj/breakup_line.o obj/strtostr.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

bin/repair_sequence2: obj/breakup_line.o obj/open_compressed.o obj/repair_sequence2.o obj/strtostr.o obj/write_fork.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

bin/repair_sequence3: obj/breakup_line.o obj/open_compressed.o obj/repair_sequence3.o obj/strtostr.o obj/write_fork.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

bin/compress_blat: obj/compress_blat.o obj/open_compressed.o obj/breakup_line.o obj/strtostr.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

bin/screen_blat: obj/screen_blat.o obj/breakup_line.o obj/strtostr.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

bin/mask_repeats_hashn: obj/breakup_line.o obj/get_name.o obj/hashn.o obj/hist_lib_hashn.o obj/mask_repeats_hashn.o obj/next_prime.o obj/open_compressed.o obj/pattern.o obj/read.o obj/read_file.o obj/strtostr.o obj/time_used.o obj/write_fork.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

bin/histogram_hashn: obj/get_name.o obj/hashn.o obj/hist_lib_hashn.o obj/histogram_hashn.o obj/next_prime.o obj/open_compressed.o obj/pattern.o obj/read.o obj/read_file.o obj/time_used.o obj/write_fork.o obj/breakup_line.o obj/strtostr.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

bin/histogram_hashl: obj/hashl.o obj/hashl_metadata.o obj/histogram_hashl.o obj/next_prime.o obj/open_compressed.o obj/time_used.o obj/write_fork.o obj/breakup_line.o obj/strtostr.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

bin/check_barcodes: obj/check_barcodes.o obj/breakup_line.o obj/strtostr.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

bin/arachne_create_xml: obj/arachne_create_xml.o obj/open_compressed.o obj/parse_readnames.o obj/write_fork.o obj/breakup_line.o obj/strtostr.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

bin/extract_seq_and_qual: obj/extract_seq_and_qual.o obj/breakup_line.o obj/open_compressed.o obj/pattern.o obj/strtostr.o obj/write_fork.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

bin/split_fasta: obj/split_fasta.o obj/open_compressed.o obj/write_fork.o obj/breakup_line.o obj/strtostr.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

bin/sort_blast: obj/sort_blast.o obj/open_compressed.o obj/breakup_line.o obj/strtostr.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

ifeq ($(OS), Linux)

bin/mask_repeats_hashz: obj/breakup_line.o obj/open_compressed.o obj/get_name.o obj/hashz.o obj/hist_lib_hashz.o obj/mask_repeats_hashz.o obj/next_prime.o obj/pattern.o obj/read.o obj/read_lib.o obj/strtostr.o obj/time_used.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS) -lgmp

bin/histogram_hashz: obj/open_compressed.o obj/get_name.o obj/hashz.o obj/hist_lib_hashz.o obj/histogram_hashz.o obj/next_prime.o obj/pattern.o obj/read.o obj/read_lib.o obj/time_used.o obj/breakup_line.o obj/strtostr.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS) -lgmp

bin/barcode_separation: obj/barcode_separation.o obj/breakup_line.o obj/open_compressed.o obj/strtostr.o obj/write_fork.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

bin/barcode_separation2: obj/barcode_separation2.o obj/breakup_line.o obj/open_compressed.o obj/strtostr.o obj/write_fork.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

depend/split_bam.d: split_bam.cc
	$(CXX) -MM $(CPPFLAGS) $(CXXFLAGS) $(PACBIO_INCLUDES) $< | sed 's,\($*\)\.o[ :]*,obj/\1.o $@ : ,g' > $@
obj/split_bam.o: split_bam.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $(PACBIO_INCLUDES) -o $@ $<
bin/split_bam: obj/split_bam.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS) $(PACBIO_LIBS)

depend/filter_bam.d: filter_bam.cc
	$(CXX) -MM $(CPPFLAGS) $(CXXFLAGS) $(PACBIO_INCLUDES) $< | sed 's,\($*\)\.o[ :]*,obj/\1.o $@ : ,g' > $@
obj/filter_bam.o: filter_bam.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $(PACBIO_INCLUDES) -o $@ $<
bin/filter_bam: obj/filter_bam.o obj/open_compressed.o obj/breakup_line.o obj/strtostr.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS) $(PACBIO_LIBS)

depend/extract_good_read_names.d: extract_good_read_names.cc
	$(CXX) -MM $(CPPFLAGS) $(CXXFLAGS) $(PACBIO_INCLUDES) $< | sed 's,\($*\)\.o[ :]*,obj/\1.o $@ : ,g' > $@
obj/extract_good_read_names.o: extract_good_read_names.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $(PACBIO_INCLUDES) -o $@ $<
bin/extract_good_read_names: obj/extract_good_read_names.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS) $(PACBIO_LIBS)

depend/filter_bam_alignments.d: filter_bam_alignments.cc
	$(CXX) -MM $(CPPFLAGS) $(CXXFLAGS) $(PACBIO_INCLUDES) $< | sed 's,\($*\)\.o[ :]*,obj/\1.o $@ : ,g' > $@
obj/filter_bam_alignments.o: filter_bam_alignments.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $(PACBIO_INCLUDES) -o $@ $<
bin/filter_bam_alignments: obj/filter_bam_alignments.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS) $(PACBIO_LIBS)

depend/pacbio_read_stats.d: pacbio_read_stats.cc
	$(CXX) -MM $(CPPFLAGS) $(CXXFLAGS) $(PACBIO_INCLUDES) $< | sed 's,\($*\)\.o[ :]*,obj/\1.o $@ : ,g' > $@
obj/pacbio_read_stats.o: pacbio_read_stats.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $(PACBIO_INCLUDES) -o $@ $<
bin/pacbio_read_stats: obj/pacbio_read_stats.o obj/breakup_line.o obj/open_compressed.o obj/pretty_print.o obj/strtostr.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS) $(PACBIO_LIBS)

depend/add_passes.d: add_passes.cc
	$(CXX) -MM $(CPPFLAGS) $(CXXFLAGS) $(PACBIO_INCLUDES) $< | sed 's,\($*\)\.o[ :]*,obj/\1.o $@ : ,g' > $@
obj/add_passes.o: add_passes.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $(PACBIO_INCLUDES) -o $@ $<
bin/add_passes: obj/add_passes.o obj/open_compressed.o obj/write_fork.o obj/breakup_line.o obj/strtostr.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS) $(PACBIO_LIBS)

depend/add_quality.d: add_quality.cc
	$(CXX) -MM $(CPPFLAGS) $(CXXFLAGS) $(PACBIO_INCLUDES) $< | sed 's,\($*\)\.o[ :]*,obj/\1.o $@ : ,g' > $@
obj/add_quality.o: add_quality.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $(PACBIO_INCLUDES) -o $@ $<
bin/add_quality: obj/add_quality.o obj/open_compressed.o obj/write_fork.o obj/breakup_line.o obj/strtostr.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS) $(PACBIO_LIBS)

depend/extract_bam_well_sizes.d: extract_bam_well_sizes.cc
	$(CXX) -MM $(CPPFLAGS) $(CXXFLAGS) $(PACBIO_INCLUDES) $< | sed 's,\($*\)\.o[ :]*,obj/\1.o $@ : ,g' > $@
obj/extract_bam_well_sizes.o: extract_bam_well_sizes.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $(PACBIO_INCLUDES) -o $@ $<
bin/extract_bam_well_sizes: obj/extract_bam_well_sizes.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS) $(PACBIO_LIBS)

depend/find_kmers.d: find_kmers.cc
	$(CXX) -MM $(CPPFLAGS) $(CXXFLAGS) $< | sed 's,\($*\)\.o[ :]*,obj/\1.o $@ : ,g' > $@
obj/find_kmers.o: find_kmers.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) -o $@ $<
bin/find_kmers: obj/find_kmers.o obj/open_compressed.o obj/hashp.o obj/next_prime.o obj/breakup_line.o obj/strtostr.o obj/time_used.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS) -lpthread

else

depend/add_passes.d:
	touch $@
depend/find_kmer.d:
	touch $@

endif

# ===========================================
# Maintenance section - dependencies, cleanup
# ===========================================

.PHONY: clean distclean depend

# only include previously made depend files; only create depend files
# if explicitly told to do so (updating of existing ones is automatic,
# though); if you want dependencies, run make depend once, then whenever
# a new source file is added

clean:
	-rm -f core *.core obj/*.o *[~%] include/*[~%]

distclean: clean
	-rm -f depend/*.d

EXISTING_DEPS := $(wildcard depend/*.d)

ifneq ($(EXISTING_DEPS),)
include $(EXISTING_DEPS)
endif

SRCS := $(wildcard *.cc)
DEPS := $(SRCS:%.cc=depend/%.d)

depend: $(DEPS)

depend/%.d: %.cc
	$(CXX) -MM $(CPPFLAGS) $(CXXFLAGS) $< | sed 's,\($*\)\.o[ :]*,obj/\1.o $@ : ,g' > $@
