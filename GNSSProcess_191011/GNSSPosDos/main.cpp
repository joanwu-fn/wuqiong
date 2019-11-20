#include "../Common/BaseFunction.h"
#include "../pNetcontrol/pNetcontrol.h"
#include <Windows.h>


int main(int argc, char *argv[])
{
//   	RePutUPD();
//     return 0;

	Netcontrol NetProcess;
	if(argc != 2)
	{
		printf("ctrl_file address is not well given\n");
		Sleep(2000);
		return 0;
	}
	else
	{
		NetProcess.ReadCtl(argv[1]);
		NetProcess.IniNet();
		char output[200];
		sprintf(output,"%s\\endoutput\\Statistic.txt",NetProcess.m_Option.Directory.ProDir);
		NetProcess.m_ftstatis = fopen(output,"w");

		int MJDN;
		for(MJDN = NetProcess.m_Option.StartTime.MJDN; MJDN <= NetProcess.m_Option.EndTime.MJDN; MJDN++)
		{

// 			NetProcess.ProcessCompareObs(MJDN);

			printf("%6d\n",MJDN);
			NetProcess.InitailFile(MJDN);			
			NetProcess.Process(MJDN);	
			NetProcess.Release();
		}

		NetProcess.OutPutUPDStatis();
		fclose(NetProcess.m_ftstatis);
		
	}
	return 1;
}