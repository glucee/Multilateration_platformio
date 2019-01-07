#include <Arduino.h>
#include <Multilateration.h>
using namespace std;

void setup() {
    Serial.begin(9600);
    vector <multilateration::Point<float>> stations = {
    multilateration::Point<float>  (1,1),
    multilateration::Point<float>  (1,0),
    multilateration::Point<float>  (0,1),
    multilateration::Point<float>  (0,0)};
    vector <float> distances = {0.1, 0.5, 0.5, 1.3};
    multilateration::Point<float> source1 = multilateration::findSourceLocation (distances, stations, multilateration::Point<float> (6.8, 6.8));
    multilateration::Point<float> source2 = multilateration::findSourceLocation (distances, stations);
    source1 = multilateration::circularClipping(source1, multilateration::Point<float> (1, 1), (float)0.3);
    Serial.print(source1.x,6);
    Serial.print (" ");
    Serial.print(source1.y,6);
    Serial.print ("     ");
    Serial.print(source2.x,6);
    Serial.print(" ");
    Serial.print(source2.y,6);
    
}
void loop() {
  // put your main code here, to run repeatedly:
  
}