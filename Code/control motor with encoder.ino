/* */

// Include the AccelStepper library:
#include <AccelStepper.h>
#include <Wire.h>
#include <AS5600.h>

#ifdef ARDUINO_SAMD_VARIANT_COMPLIANCE
  #define Serial SerialUSB
  #define SYS_VOL   3.3
#else
  #define Serial Serial
  #define SYS_VOL   5
#endif

// Define stepper motor connections and motor interface type. Motor interface type must be set to 1 when using a driver:
#define dirPin 2
#define stepPin 3
#define enablePin 4
#define motorInterfaceType 1


// Create a new instance of the AccelStepper class:
AccelStepper stepper = AccelStepper(motorInterfaceType, stepPin, dirPin);
AMS_5600 ams5600;



float zero = 82.2149960000-0.9570007300-0.0870010400



;

float desired = 0; // Default needs to change to the perpendicular position not yet measured

float maxAngle = 30;
float minAngle = -10;
float anglePres = 0.07;


float angle;
float prevAngle;

float convertRawAngleToDegrees(word newAngle)
{
  /* Raw data reports 0 - 4095 segments, which is 0.087 of a degree */
  float retVal = newAngle * 0.087;
  return retVal;
}


void setup() {
  Serial.begin(9600); // Start Serial Monitor at Board Rate 9600
  Serial.println("Setup Initiated");

  stepper.setEnablePin(enablePin);
  stepper.setPinsInverted(false, false, true);

  stepper.enableOutputs();        //Enables Motor Pin Outputs
  stepper.setMaxSpeed(1000);
  stepper.setAcceleration(5000);
  stepper.disableOutputs();
 
  Wire.begin();
  if(ams5600.detectMagnet() == 0 ){
    while(1){
        if(ams5600.detectMagnet() == 1 ){
            Serial.print("Current Magnitude: ");
            Serial.println(ams5600.getMagnitude());
            break;
        }
        else{
            Serial.println("Can not detect magnet");
        }
        delay(1000);
    }
  }
  ams5600.setStartPosition(minAngle);
  ams5600.setEndPosition(maxAngle);  
  Serial.println("Setup Complete");
} 

void loop() {
  if (Serial.available() > 0){
    
    float desireduncheck = Serial.parseFloat(SKIP_ALL);
    
    if (desireduncheck < maxAngle && desireduncheck > minAngle && desireduncheck != 0){
      desired = desireduncheck;
      Serial.println("New angle: "+String( desired,DEC));
    }
    else{Serial.println("Out of range");} 
  }

  angle = convertRawAngleToDegrees(ams5600.getRawAngle())-zero;

  
  if (abs(angle-desired) > anglePres){

    stepper.enableOutputs();              //Enables Motor Pin Outputs
    
    if (angle > desired){
      
      stepper.setSpeed(-500);
      stepper.runSpeed(); 
    }
    else if (angle < desired){
      stepper.setSpeed(500);
      stepper.runSpeed();
    }
    Serial.println(String(angle,DEC));
    
  }
  else if (abs(angle-prevAngle) > anglePres){
    Serial.println(String(angle,DEC));  
  }
  
  stepper.disableOutputs();             //Disables Motor Pin Outputs
  prevAngle = angle;
}
