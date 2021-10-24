#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>


/*                                          STRUCT                                                         */
struct Data {
	/*struct to maintain the value as well an index of specific dataset*/
	int index;
	float dataValue;
};
/*                                     PRINTING FUNCTIONS                                                   */
void print(struct Data** array, int size) {
	/*print function to print the required values of the interpolation table*/
	int i = 0;
	for (; i < size; i++) {
		int j = 0;
		for (; j < size + 1; j++) {
			if (array[i][j].dataValue != INT_MIN)       //checking if the current is one of the used indexes
				printf("%15f", array[i][j].dataValue);
			else
				printf("   ");
		}
		printf("\n");
	}
}
void printLine(int n) {
	/*helper function to print line for interface*/
	int i = 0;
	for (; i < n; i++)
		printf("-");
	printf("\n");
}
/*                                  HELPER FUNCTIONS                                                         */
void setValues(int index, float dataValue, struct Data* arr) {
	/*helper function to set the values of a particual object*/
	arr->index = index;
	arr->dataValue = dataValue;
}
void setArrayTo0(struct Data** array, int size) {
	/*helper function to initialize the Data array*/
	int i = 0;
	for (; i < size; i++) {
		int j = 0;
		for (; j < size + 1; j++)
			array[i][j].dataValue = INT_MIN;  //storing INT_MIN at each index
	}
}
bool equalOrUnequal(struct Data** arr, int size) {
	/*helper function to check whether the given dataset has equal intervals or unequal intervals*/
	bool equal = true;
	int xValue1 = 0;
	int xValue2 = 1;
	float difference = arr[xValue2][0].dataValue - arr[xValue1][0].dataValue;
	while (xValue2 < size) {
		if (difference != arr[xValue2][0].dataValue - arr[xValue1][0].dataValue) {
			//if even once the difference comes out ot be different, the dataset is considered unequal
			equal = false;
			break;
		}
		xValue1++;
		xValue2++;
	}
	return equal;
}
struct Data* colToRow(struct Data** arr, int col, int row) {
	/*helper function to convert a col into a 1d array*/
	struct Data* array = (struct Data*)malloc(sizeof(struct Data) * (row));
	int i = 0;
	for (; i < row; i++)
		array[i] = arr[i][col];
	return array;
}
void subtractEqual(struct Data* arr, struct Data** array, int col, int size) {
	/*helper function to calculate a specific coloum*/
	int yValue1 = 0;
	int yValue2 = 1;
	while (yValue2 <= size) {
		//a specific coloum entity is calculated by subtracting corresponding entities of the previous coloum
		array[yValue1][col].dataValue = arr[yValue2].dataValue - arr[yValue1].dataValue;
		array[yValue1][col].index = yValue1;
		yValue1++;
		yValue2++;
	}
}
void subtractUnequal(struct Data* arr, struct Data** array, int col, int size, int xValue1, int xValue2) {
	/*This function is used by newton divided difference formula*/
	int yValue1 = 0;
	int yValue2 = 1;
	while (yValue2 <= size) {
		/*a specific coloum entity is calculated by subtracting corresponding
		entities of the previous coloum and dividing them by the range they lie in*/
		array[yValue1][col].dataValue = (arr[yValue2].dataValue - arr[yValue1].dataValue) /
			(array[xValue2][0].dataValue - array[xValue1][0].dataValue);
		array[yValue1][col].index = yValue1;
		yValue1++;
		yValue2++;
		xValue2++;
		xValue1++;
	}
}
void createTableEqual(struct Data** arr, int rows, int col) {
	/*helper function to set the values of the table*/
	int i = 1;
	for (; i < col; i++) {
		struct Data* array = (struct Data*)malloc(sizeof(struct Data) * (rows));
		array = colToRow(arr, i, rows);             //convert col to row
		subtractEqual(array, arr, i + 1, rows - i); //set the values by taking differences
	}
}
void createTableUnequal(struct Data** arr, int rows, int col) {
	/*helper function to create table for unequal difference*/
	int i = 1;
	int xValue1 = 0;
	int xValue2 = 1;
	for (; i < col; i++) {
		struct Data* array = (struct Data*)malloc(sizeof(struct Data) * (rows));
		array = colToRow(arr, i, rows);  //converting col to 1d array
		subtractUnequal(array, arr, i + 1, rows - i, xValue1, xValue2);
		xValue2++;
	}
}
float* generatePatternForP(float p, int size, bool check) {
	/*helper function to generate the P values pattern for newton forward and backward formulas*/
	int i = 1;
	float current = 0;
	float* pValueArray = (float*)malloc(sizeof(float) * size);
	int index = 0;
	for (; i <= size; i++) {
		if (check == false)   //if check == false, then create according to newton forward
			current = p - i;
		else                  //else create according to newton backward
			current = p + i;
		if (index - 1 >= 0)
			current = current * pValueArray[index - 1];
		else
			current = current * p;
		pValueArray[index] = current; //store the current value in 1d array
		index++;
	}
	return pValueArray;
}
float* generatePatternForUnequal(float interpolateValue, int size, struct Data** arr, int exclude) {
	/*helper function to generate the P values pattern for newton forward and backward formulas*/
	int i = 0;
	float* xValueArray = (float*)malloc(sizeof(float) * size);
	int index = 0;
	for (; i < size; i++) {
		if (i != exclude) { //exclude the specific term that is not part of the formula
			float current = interpolateValue - arr[i][0].dataValue;
			if (index - 1 >= 0)
				current = current * xValueArray[index - 1];
			xValueArray[index] = current;
			index++;
		}
	}
	return xValueArray;
}
float* generatePatternForGF(float p, int size, bool check) {
	/*helper function to generate pattern for gauss forward and backward formulas*/
	int i = 0;
	bool alternate = true;
	int minus = 1;
	int plus = 1;
	float* pValueArray = (float*)malloc(sizeof(float) * size);
	int index = 0;
	float current = 0;
	for (; i < size; i++) {
		if (alternate == check) {
			current = p - minus;
			minus++;
			if (check == false)  //if check == false calculate for gauss backward else for gauss forward
				alternate = true;
			else
				alternate = false;
		}
		else {
			current = p + plus;
			plus++;
			if (check == false)
				alternate = false;
			else
				alternate = true;
		}
		if (index - 1 >= 0)
			current = current * pValueArray[index - 1];
		else
			current = current * p;

		pValueArray[index] = current;  //store the current answer in the array
		index++;
	}
	return pValueArray;
}
int* factorial(int size) {
	/*helper function to generate factorial*/
	int i = 3;
	int* array = (int*)malloc(sizeof(int) * size);
	array[0] = 2;
	int index = 1;
	for (; index < size; index++) {
		array[index] = array[index - 1] * i; //store the results in 1d array
		i++;
	}
	return array;
}
int findX0Index(struct Data** array, int row, float interpolateValue, float lowerlimit, float upperlimit, float* p) {
	/*helper function to calculate the x0 based on p value*/
	float diff = array[1][0].dataValue - array[0][0].dataValue; //calculate 'h'
	int i = 0;
	for (; i < row; i++) {
		*p = (interpolateValue - array[i][0].dataValue) / diff;
		if (*p >= lowerlimit && *p <= upperlimit) //if p lies in limit then return current index
			return i;
	}
	return 0;
}
void setIndex(struct Data** array, int smallestPValue, int row, int col) {
	/*helper function to set the indexes according to the x0*/
	int j = 0;
	for (; j < col; j++) {
		int i = 0;
		int currentIndex = smallestPValue;
		for (; i < row; i++) {
			array[i][j].index = currentIndex; //set index
			currentIndex++;
		}
	}
}
float returnYValue(struct Data** arr, int row, int col, int arrIndex, bool* breakLoop) {
	/*helper function to return the value of a specific index based on the objects index attribute*/
	int i = 0;
	for (; i < row - col + 1; i++) {
		if (arr[i][col].index == arrIndex) {
			*breakLoop = false;
			return arr[i][col].dataValue;
		}
	}
	*breakLoop = true;
	return -1;
}
/*                                      EQUAL FORMULAS                                                      */
float newtonForward(struct Data** array, int row, int col, float interpolateValue) {
	/*function to interpolate using newton forward difference formula*/
	float diff = array[1][0].dataValue - array[0][0].dataValue;  //calculate 'h'
	float p = (interpolateValue - array[0][0].dataValue) / diff; //calculate 'p'
	float* arr = generatePatternForP(p, row - 2, false);         //generate p pattern
	int* farr = factorial(row - 2);                              //generate factorials
	float answer = array[0][1].dataValue + p * array[0][2].dataValue;
	int i = 3;
	int index = 0;
	for (; i < col; i++) {
		answer += (arr[index] * array[0][i].dataValue) / farr[index]; //update the answer on each calculation
		index++;
	}
	return answer;
}
float newtonBackward(struct Data** array, int row, int col, float interpolateValue) {
	/*function to interpolate using newton backward difference formula*/
	float diff = array[1][0].dataValue - array[0][0].dataValue;        //calculate 'h'
	float p = (interpolateValue - array[row - 1][0].dataValue) / diff; //calculate 'p'
	float* arr = generatePatternForP(p, row - 2, true);                //generate p pattern
	int* farr = factorial(row - 2);                                    //generate factorials
	float answer = array[row - 1][1].dataValue + p * array[row - 2][2].dataValue;
	int i = 3;
	int index = 0;
	for (; i < col; i++) {
		answer += (arr[index] * array[row - i][i].dataValue) / farr[index]; //update the answer on each calculation
		index++;
	}
	return answer;
}
float gaussForward(struct Data** array, int row, int col, float interpolateValue) {
	/*function to interpolate using gauss forward formula*/
	float p = 0;
	bool breakLoop = false;
	int cindex = -1;
	int count = 0;
	int index = findX0Index(array, row, interpolateValue, 0, 1, &p); //finding the x0
	setIndex(array, -1 * index, row, col);                           //setting index based on x0
	float* pArr = generatePatternForGF(p, row - 2, true);            //generate pattern
	int* farr = factorial(row - 2);                                  //generate factorials
	float answer = returnYValue(array, row, 1, 0, &breakLoop);
	if (breakLoop == true)
		return 0;
	float temp = p * returnYValue(array, row, 2, 0, &breakLoop);
	if (breakLoop == true)
		return answer;
	answer += temp;
	int i = 3;
	index = 0;
	for (; i < col && breakLoop == false; i++) {
		temp = returnYValue(array, row, i, cindex, &breakLoop);
		count++;
		if (count % 2 == 0)
			cindex = cindex - 1;
		if (breakLoop == true)  //if the required index does not exist, return from function 
			return answer;
		answer += temp * pArr[index] / farr[index];  //update answer on each calculation
		index++;
	}
	return answer;
}
float gaussBackward(struct Data** array, int row, int col, float interpolateValue) {
	/*function to interpolate using gauss backward formula*/
	float p = 0;
	bool breakLoop = false;
	int cindex = -1;
	int count = 1;
	int index = findX0Index(array, row, interpolateValue, -1, 0, &p); //find the x0 index
	setIndex(array, -1 * index, row, col);                            //set the indexes based on x0
	float* pArr = generatePatternForGF(p, row - 2, false);            //generate pattern for p
	int* farr = factorial(row - 2);                                   //generate factorials
	float answer = returnYValue(array, row, 1, 0, &breakLoop);
	if (breakLoop == true)
		return 0;
	float temp = p * returnYValue(array, row, 2, -1, &breakLoop);
	if (breakLoop == true)
		return answer;
	answer += temp;
	int i = 3;
	index = 0;
	for (; i < col && breakLoop == false; i++) {
		temp = returnYValue(array, row, i, cindex, &breakLoop);
		count++;
		if (count % 2 == 0)
			cindex = cindex - 1;
		if (breakLoop == true)  //if the required index does not exist, return from function
			return answer;
		answer += temp * pArr[index] / farr[index]; //update answer on each calculation
		index++;
	}
	return answer;
}
float sterling(struct Data** array, int row, int col, float interpolateValue) {
}
/*                                      UNEQUAL FORMULAS                                                    */
float lagrange(struct Data** array, int row, float interpolateValue) {
	/*function to interpolate using lagrange formula*/
	float answer = 0;
	int i = 0;
	for (; i < row; i++) {
		float* numerator = generatePatternForUnequal(interpolateValue, row, array, i); //calculate numerator
		float* denominator = generatePatternForUnequal(array[i][0].dataValue, row, array, i); //calcultae denominator
		answer += ((numerator[row - 2] / denominator[row - 2]) * array[i][1].dataValue); //update answer on each calculation
	}
	return answer;
}
float nfDividedDifference(struct Data** array, int row, int col, float interpolateValue) {
	/*function to interpolate using newton divided difference formula*/
	float* arr = generatePatternForUnequal(interpolateValue, row, array, -1); //generate pattern for p
	float answer = array[0][1].dataValue;
	int i = 2;
	int index = 0;
	for (; i < col; i++) {
		answer += (arr[index] * array[0][i].dataValue);  //update answer on each calculation
		index++;
	}
	return answer;
}
int main() {
	system("clear");
	int totalDataSet;
	printf("%s", "Enter the number of dataset values: ");
	scanf("%d", &totalDataSet);  //get input for total number of datasets
	struct Data* obj[totalDataSet];
	int i = 0;
	for (; i < totalDataSet; i++) {
		obj[i] = (struct Data*)malloc(sizeof(struct Data) * (totalDataSet + 1)); //allocate memory
	}
	setArrayTo0(obj, totalDataSet); //initialize array
	i = 0;
	for (; i < totalDataSet; i++) {
		float xValue;
		float yValue;
		printf("%s", "Enter value of x: ");
		scanf("%f", &xValue);
		printf("%s", "Enter value of y: ");
		scanf("%f", &yValue);
		setValues(i, xValue, &obj[i][0]);   //set x value
		setValues(i, yValue, &obj[i][1]);   //set y value
	}
	float interpolateValue;
	printf("%s", "Enter the value on which you want to interpolate: ");
	scanf("%f", &interpolateValue);   //take input  for interpolation value
	print(obj, totalDataSet);
	bool h = equalOrUnequal(obj, totalDataSet);    //check to see if dataset is unequal or equal
	if (h == true) { //if bool == true solve for equal dataset
		system("clear");
		printLine(100);
		createTableEqual(obj, totalDataSet, totalDataSet + 1); //create table for equal dataset
		print(obj, totalDataSet);
		printLine(100);
		printf("\n\n");
		printf("%s", "Equal Interval Formulas\n");
		printLine(50);
		float nfAnswer = newtonForward(obj, totalDataSet, totalDataSet + 1, interpolateValue); //newton forward
		printf("%s %f\n\n", "Newton Forward Difference Formula:  ", nfAnswer);
		printLine(50);
		float nbAnswer = newtonBackward(obj, totalDataSet, totalDataSet + 1, interpolateValue); //newton backward
		printf("%s %f\n\n", "Newton Backward Difference Formula: ", nbAnswer);
		printLine(50);
		float lagrangeAnswer = lagrange(obj, totalDataSet, interpolateValue);                   //lagrange formula
		printf("%s %f\n\n", "Lagrange Formula:                   ", lagrangeAnswer);
		printLine(50);
		float gfAnswer = gaussForward(obj, totalDataSet, totalDataSet + 1, interpolateValue);   //gauss forward
		printf("%s %f\n\n", "Gauss Forward Difference Formula:   ", gfAnswer);
		printLine(50);
		float gbAnswer = gaussBackward(obj, totalDataSet, totalDataSet + 1, interpolateValue);  //gauss backward
		printf("%s %f\n\n", "Gauss Backward Difference Formula:  ", gbAnswer);
		printLine(50);
	}
	else {
		system("clear");
		printLine(100);
		createTableUnequal(obj, totalDataSet, totalDataSet + 1);  //create table for unequal dataset
		print(obj, totalDataSet);
		printLine(100);
		printf("\n\n");
		printf("%s", "Unequal Interval Formulas\n");                
		printLine(50);
		float nfddAnswer = nfDividedDifference(obj, totalDataSet, totalDataSet + 1, interpolateValue); //newton divided difference
		printf("%s %f\n\n", "Newton Divided Difference Formula: ", nfddAnswer);
		printLine(50);
		float lagrangeAnswer = lagrange(obj, totalDataSet, interpolateValue);
		printf("%s %f\n\n", "Lagrange Formula:                  ", lagrangeAnswer);                    //lagrange formula
		printLine(50);
	}
	return 0;
}