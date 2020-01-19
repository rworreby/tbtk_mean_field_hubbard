echo "Running sim for 10x10 ..."
make gnfh > test.txt
for i in {1..10}
do
	time make gnfh > test.txt
done
echo "Running sim for 11x11 ..."
make gnfh11 > test.txt
for i in {1..10}
do
	time make gnfh11 > test.txt
done
echo "Running sim for 12x12 ..."
make gnfh12 > test.txt
for i in {1..10}
do
	time make gnfh12 > test.txt
done
