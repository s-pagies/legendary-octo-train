	    float gx_2 = gridpoints_x*0.5
            if(diffusion==true)
            {
	      if(x_pos[n]<=2)
	      {
		if(y_pos[n]>=a && y_pos<=gridpoints_y-a)
		{
		  y_pos[n] += -v;
		}
	      }
	      if(x_pos[n]>2 && x_pos[n]<gx_2-2)
	      {
		if(y_pos[n]>=a-1 && y_pos[n]<=a+1)
		{
		  x_pos[n] += v;
		}
		if(y_pos[n]>=gridpoints_y-a-1 && y_pos[n]<=gridpoints_y-a+1)
		{
		  x_pos[n] += v;
		}
	      }
	      if(x_pos[n]>=gx_2-2 && x_pos[n]<=gx_2)
	      {
		if(y_pos[n]>=a && y_pos[n]<=gridpoints_y-a)
		{
		  y_pos[n] += v;
		}
	      }
	      if(x_pos[n]>gx_2 && x_pos[n]<= gx_2+2)
	      {
		if(y_pos[n]>=a && y-pos[n]<=gridpoints_y-a)
		{
		  y_pos[n] += v;
		}
	      }
	      if(x_pos[n]>gx_2+2 && x_pos[n]<gridpoints_x-2)
	      {
		if(y_pos[n]>=a-1 && y_pos[n]<=a+1)
		{
		  x_pos[n] += -v;
		}
		if(y_pos[n]>=gridpoints_y-a-1 && y_point[n]<=gridpoints_y-a+1)
		{
		  x_pos[n] += v;
		}
	      }
	      if(x_pos[n]>=gridpoints_x-2)
	      {
		if(y_pos[n]>=a && y_pos[n]<=gridpoints_y-a)
		{
