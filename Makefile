imshrinker:
	cd src
	g++ *.cpp -o imshrinker
debug:
	g++ -std=c++20 src/Application.cc src/BitIOStream.cc src/Decoder.cc src/Encoder.cc src/Exception.cc src/FileIOStream.cc src/Image.cc src/main.cc src/SPIHT.cc -DDEBUG -Wall -g -o imshrinker
clean:
	rm -f imshrinker
