test_methods: test_methods.cpp qr.cpp
	g++ -o test_methods test_methods.cpp qr.cpp -std=c++11 -O3 -larmadillo

timed_test: timed_test.cpp qr.cpp
	g++ -o timed_test timed_test.cpp qr.cpp utils/timer.cpp -std=c++11 -O3 -larmadillo

sparse_test: sparse_test.cpp qr.cpp
	g++ -std=c++11 -O3 -o sparse_test sparse_test.cpp qr.cpp -larmadillo

clean:
	rm -f test_methods timed_test sparse_test