runTime.setDeltaT( forceTimes[(runTime.timeIndex())%numSteps+1] - forceTimes[(runTime.timeIndex())%numSteps] );

// Info<< "\t\t Dt = " << runTime.timeIndex() << endl; 

// Info<< "\t\t Dt = " << forceTimes[(runTime.timeIndex()-1)%numSteps] << endl; 