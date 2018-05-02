mkdir -p bin
for file in *.cpp; do
	g++ $file -o bin/${file%.*}
done
(cd ACO && make)
mv ACO/main bin/mmas
