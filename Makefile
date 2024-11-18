CC = g++

CPPFLAGS += -std=c++20 -Wall -g


test_vectors:
	rm -f obj/*.obj
	rm -f bin/tests
	$(CC) $(CPPFLAGS) tests/test_vectors.cpp -o bin/test_vectors
	./bin/test_vectors

test_matrices:
	rm -f obj/*.obj
	rm -f bin/tests
	$(CC) $(CPPFLAGS) tests/test_matrices.cpp -o bin/test_matrices
	./bin/test_matrices


test_LA:
	rm -f bin/tests
	$(CC) $(CPPFLAGS) tests/test_linalg.cpp -o bin/testLA
	./bin/testLA

test_LR:
	rm -f bin/tests
	$(CC) $(CPPFLAGS) tests/test_linear_regression.cpp -o bin/testLR
	./bin/testLR