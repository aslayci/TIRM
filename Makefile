all:
	g++ src_TIRM/MinRegretFramework.cc -O3 src_TIRM/utils.cc src_TIRM/anyoption.cc src_TIRM/sfmt/SFMT.c src_TIRM/MinRegretAllocator.cc src_TIRM/TimGraph.cc -o TIRM
