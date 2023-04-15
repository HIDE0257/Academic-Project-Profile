module Final_Project (input accel, speed, clk, reset, distance, brake, swcruise, output reg [1:0] led, output reg cruise, output reg warning, output reg [6:0] seg1, seg2, seg3, seg4); 
							//  button1  button2     switch0 switch9  switch1 switch2     		2-bit-led		1-bit-led called cruise_control	//1-bit-led		seg1&seg2 for current speed seg3&seg4 for desire speed
	
	reg new_clk;
	parameter max_count = 10000000;
	reg [24:0] count;
	reg [6:0] present_state, next_state, speedinput, speedout;
	reg prep = 1;

//Clock
	always @ (posedge clk) begin 
		if (count < max_count) begin 
			count <= count + 1;
		end else begin 
			count <= 0;
			new_clk <= ~new_clk;
		end
	end
	
// FF (D-FF)
	always @ (posedge new_clk) begin 
		if (reset == 1) begin //when the reset switch is on, present_state and speedinput are going to be 0;
			present_state = 0; //it is for a current speed
			speedinput = 0; // it is for a desired speed
		end else begin 
			present_state <= next_state; //when the reset switch is off, both of them are variable 
			speedinput <= speedout;
		end 
	end 
	
	//reset --> switch 0
	//speed --> switch 1
	//cruise --> switch 2
	//distance --> switch 9
	
	//accel --> button 1
	//brake --> button 2
	
	
// next_state
	always @ (*) begin 
		if (distance == 0) begin //distance (sensor) is the highest priority; when it's off, all the orders below should be implemented
			if (brake == 1) begin
				if (swcruise == 0 && prep == 1) begin  //when the cruise control is off, the orders below  should be implemented
					cruise = 1;
					if (speed == 1) begin //when the speed switch is on, 
						if (speedinput <= 64) begin
							speedout = speedinput + 1; //speedout is going to increase by 1 until speedinput reaches at 64
						end else begin 	// when speedinput is over 65,
							speedout = speedinput;	//speedout stays the same as speedinput
						end
					end else begin		//when the speed switch is off,
						speedout = speedinput;	//speedout stays the same as speedinput
					end
				 
					if (speed == 0 && accel == 1) begin 	//When both the speed switch and the gas pedal are off, 
						if (present_state < speedinput) begin //
							next_state = present_state + 1;
						end else if (present_state > speedinput) begin 
							next_state = present_state - 1;
						end else begin 
							next_state = present_state;
						end
					end else if (speed == 0 && accel == 0) begin 
						if (present_state <= 64) begin 
							next_state = present_state + 1;
						end else begin 
							next_state = present_state;
						end 
					end 
				end else if(prep == 1 && swcruise == 1) begin
					cruise = 0;
					if (accel == 0) begin
						if (present_state <= 64) begin 
							next_state = present_state + 1;	
						end else begin
							next_state = present_state;
						end
					end else if (brake == 0) begin 
						if (present_state > 0) begin 
							next_state = present_state - 1;
						end else begin 
							next_state = present_state;
						end
					end
				
				end else if (prep == 0 && swcruise == 0) begin
					cruise = 0;
					if (accel == 0) begin
						if (present_state <= 64) begin 
							next_state = present_state + 1;	
						end else begin
							next_state = present_state;
						end
					end else if (brake == 0) begin 
						if (present_state > 0) begin 
							next_state = present_state - 1;
						end else begin 
							next_state = present_state;
						end
					end
				
				end else begin
					prep = 1;
					cruise = 1;
				end
				
			end else begin
				prep = 0;	
				if (present_state > 0) begin 
					next_state = present_state - 1;
				end else begin 
					next_state = present_state;
				end
			
				if (accel == 0) begin
					if (present_state <= 64) begin 
						next_state = present_state + 1;	
					end else begin
						next_state = present_state;
					end
				end
			end
			
		end else begin 
			prep = 0;
			if (present_state > 0) begin 
				next_state = present_state - 1;
			end else begin 
				next_state = present_state;
			end
		end
	end


// output
	always @(*) begin 
		if (distance == 0) begin	//when the destance switch is off, then the order below is implemented 
			warning = 1'b0;
			if (0 <= present_state && present_state < speedout) begin	//0 <= current speed < desired speed
				led = 2'b01; //tell to accelerate; the the most right LED turns on
				
			end else if (present_state == speedout) begin
				led = 2'b00; //when the current speed is the same as the desired speed, no LEDs are on
			
			end else	if (speedout < present_state && present_state <= 64) begin  //desired speed < current speed <= 64
				led = 2'b10; //tell to decelerate
		
			end 
		end else begin //when the distance switch is on, notify the cruise is decelerating
			led = 2'b10;
			warning = 1'b1;
		end
		
		
		
		
		//Display the current speed from 0 to 64.
		case (present_state)  //The display depends on the change in present_state
			00: begin
				seg1 = 7'b1000000;
				seg2 = 7'b1000000;
			end 
			
			01: begin 
				seg1 = 7'b1111001;
				seg2 = 7'b1000000;
			end
			02: begin 
				seg1 = 7'b0100100;
				seg2 = 7'b1000000;
			end
			03: begin 
				seg1 = 7'b0110000;
				seg2 = 7'b1000000;
			end 
			04: begin 
				seg1 = 7'b0011001;
				seg2 = 7'b1000000;
			end 
			05: begin 
				seg1 = 7'b0010010;
				seg2 = 7'b1000000;
			end 
			06: begin 
				seg1 = 7'b0000010;
				seg2 = 7'b1000000;
			end 
			07: begin 
				seg1 = 7'b1111000;
				seg2 = 7'b1000000;
			end
			08: begin 
				seg1 = 7'b0000000;
				seg2 = 7'b1000000;
			end 
			09: begin 
				seg1 = 7'b0010000;
				seg2 = 7'b1000000;
			end 
			10: begin
				seg1 = 7'b1000000;
				seg2 = 7'b1111001;
			end
			11: begin 
				seg1 = 7'b1111001;
				seg2 = 7'b1111001;
			end 
			12: begin 
				seg1 = 7'b0100100;
				seg2 = 7'b1111001;
			end 
			13: begin 
				seg1 = 7'b0110000;
				seg2 = 7'b1111001;
			end 
			14: begin 
				seg1 = 7'b0011001;
				seg2 = 7'b1111001;
			end 
			15: begin
				seg1 = 7'b0010010;
				seg2 = 7'b1111001;
			end 
			16: begin 
				seg1 = 7'b0000010;
				seg2 = 7'b1111001;
			end 
			17: begin 
				seg1 = 7'b1111000;
				seg2 = 7'b1111001;
			end 
			18: begin 
				seg1 = 7'b0000000;
				seg2 = 7'b1111001;
			end 
			19: begin 
				seg1 = 7'b0010000;
				seg2 = 7'b1111001;
			end 
			20: begin 
				seg1 = 7'b1000000;
				seg2 = 7'b0100100;
			end 
			21: begin 
				seg1 = 7'b1111001;
				seg2 = 7'b0100100;
			end 
			22: begin 
				seg1 = 7'b0100100;
				seg2 = 7'b0100100;
			end 
			23: begin 
				seg1 = 7'b0110000;
				seg2 = 7'b0100100;
			end 
			24: begin 
				seg1 = 7'b0011001;
				seg2 = 7'b0100100;
			end 
			25: begin
				seg1 = 7'b0010010;
				seg2 = 7'b0100100;
			end 
			26: begin 
				seg1 = 7'b0000010;
				seg2 = 7'b0100100;
			end 
			27: begin 
				seg1 = 7'b1111000;
				seg2 = 7'b0100100;
			end 
			28: begin 
				seg1 = 7'b0000000;
				seg2 = 7'b0100100;
			end 
			29: begin 
				seg1 = 7'b0010000;
				seg2 = 7'b0100100;
			end 
			30: begin 
				seg1 = 7'b1000000;
				seg2 = 7'b0110000;
			end
			31: begin 
				seg1 = 7'b1111001;
				seg2 = 7'b0110000;
			end 
			32: begin 
				seg1 = 7'b0100100;
				seg2 = 7'b0110000;
			end 
			33: begin 
				seg1 = 7'b0110000;
				seg2 = 7'b0110000;
			end 
			34: begin 
				seg1 = 7'b0011001;
				seg2 = 7'b0110000;
			end 
			35: begin
				seg1 = 7'b0010010;
				seg2 = 7'b0110000;
			end 
			36: begin 
				seg1 = 7'b0000010;
				seg2 = 7'b0110000;
			end 
			37: begin 
				seg1 = 7'b1111000;
				seg2 = 7'b0110000;
			end 
			38: begin 
				seg1 = 7'b0000000;
				seg2 = 7'b0110000;
			end 
			39: begin 
				seg1 = 7'b0010000;
				seg2 = 7'b0110000;
			end 
			40: begin 
				seg1 = 7'b1000000;
				seg2 = 7'b0011001;
			end
			41: begin 
				seg1 = 7'b1111001;
				seg2 = 7'b0011001;
			end
			42: begin 
				seg1 = 7'b0100100;
				seg2 = 7'b0011001;
			end
			43: begin 
				seg1 = 7'b0110000;
				seg2 = 7'b0011001;
			end
			44: begin 
				seg1 = 7'b0011001;
				seg2 = 7'b0011001;
			end 
			45: begin 
				seg1 = 7'b0010010;
				seg2 = 7'b0011001;
			end 
			46: begin 
				seg1 = 7'b0000010;
				seg2 = 7'b0011001;
			end
			47: begin 
				seg1 = 7'b1111000;
				seg2 = 7'b0011001;
			end 
			48: begin 
				seg1 = 7'b0000000;
				seg2 = 7'b0011001;
			end 
			49: begin 
				seg1 = 7'b0010000;
				seg2 = 7'b0011001;
			end 
			50: begin 
				seg1 = 7'b1000000;
				seg2 = 7'b0010010;
			end
			51: begin 
				seg1 = 7'b1111001;
				seg2 = 7'b0010010;
			end
			52: begin 
				seg1 = 7'b0100100;
				seg2 = 7'b0010010;
			end
			53: begin 
				seg1 = 7'b0110000;
				seg2 = 7'b0010010;
			end
			54: begin 
				seg1 = 7'b0011001;
				seg2 = 7'b0010010;
			end 
			55: begin 
				seg1 = 7'b0010010;
				seg2 = 7'b0010010;
			end 
			56: begin 
				seg1 = 7'b0000010;
				seg2 = 7'b0010010;
			end
			57: begin 
				seg1 = 7'b1111000;
				seg2 = 7'b0010010;
			end 
			58: begin 
				seg1 = 7'b0000000;
				seg2 = 7'b0010010;
			end 
			59: begin 
				seg1 = 7'b0010000;
				seg2 = 7'b0010010;
			end 
			60: begin 
				seg1 = 7'b1000000;
				seg2 = 7'b0000010;
			end
			61: begin 
				seg1 = 7'b1111001;
				seg2 = 7'b0000010;
			end
			62: begin 
				seg1 = 7'b0100100;
				seg2 = 7'b0000010;
			end
			63: begin 
				seg1 = 7'b0110000;
				seg2 = 7'b0000010;
			end
			64: begin 
				seg1 = 7'b0011001;
				seg2 = 7'b0000010;
			end
			
		endcase
		
		//Display the desired speed the user set up from 0 to 64
		case (speedinput) //The display depends on the change in speedinput
			00: begin
				seg3 = 7'b1000000;
				seg4 = 7'b1000000;
			end 
			
			01: begin 
				seg3 = 7'b1111001;
				seg4 = 7'b1000000;
			end
			02: begin 
				seg3 = 7'b0100100;
				seg4 = 7'b1000000;
			end
			03: begin 
				seg3 = 7'b0110000;
				seg4 = 7'b1000000;
			end 
			04: begin 
				seg3 = 7'b0011001;
				seg4 = 7'b1000000;
			end 
			05: begin 
				seg3 = 7'b0010010;
				seg4 = 7'b1000000;
			end 
			06: begin 
				seg3 = 7'b0000010;
				seg4 = 7'b1000000;
			end 
			07: begin 
				seg3 = 7'b1111000;
				seg4 = 7'b1000000;
			end
			08: begin 
				seg3 = 7'b0000000;
				seg4 = 7'b1000000;
			end 
			09: begin 
				seg3 = 7'b0010000;
				seg4 = 7'b1000000;
			end 
			10: begin
				seg3 = 7'b1000000;
				seg4 = 7'b1111001;
			end
			11: begin 
				seg3 = 7'b1111001;
				seg4 = 7'b1111001;
			end 
			12: begin 
				seg3 = 7'b0100100;
				seg4 = 7'b1111001;
			end 
			13: begin 
				seg3 = 7'b0110000;
				seg4 = 7'b1111001;
			end 
			14: begin 
				seg3 = 7'b0011001;
				seg4 = 7'b1111001;
			end 
			15: begin
				seg3 = 7'b0010010;
				seg4 = 7'b1111001;
			end
			16: begin 
				seg3 = 7'b0000010;
				seg4 = 7'b1111001;
			end 
			17: begin 
				seg3 = 7'b1111000;
				seg4 = 7'b1111001;
			end 
			18: begin 
				seg3 = 7'b0000000;
				seg4 = 7'b1111001;
			end 
			19: begin 
				seg3 = 7'b0010000;
				seg4 = 7'b1111001;
			end
			20: begin 
				seg3 = 7'b1000000;
				seg4 = 7'b0100100;
			end 
			21: begin 
				seg3 = 7'b1111001;
				seg4 = 7'b0100100;
			end 
			22: begin 
				seg3 = 7'b0100100;
				seg4 = 7'b0100100;
			end 
			23: begin 
				seg3 = 7'b0110000;
				seg4 = 7'b0100100;
			end 
			24: begin 
				seg3 = 7'b0011001;
				seg4 = 7'b0100100;
			end 
			25: begin
				seg3 = 7'b0010010;
				seg4 = 7'b0100100;
			end 
			26: begin 
				seg3 = 7'b0000010;
				seg4 = 7'b0100100;
			end
			27: begin 
				seg3 = 7'b1111000;
				seg4 = 7'b0100100;
			end 
			28: begin 
				seg3 = 7'b0000000;
				seg4 = 7'b0100100;
			end 
			29: begin 
				seg3 = 7'b0010000;
				seg4 = 7'b0100100;
			end 
			30: begin 
				seg3 = 7'b1000000;
				seg4 = 7'b0110000;
			end
			31: begin 
				seg3 = 7'b1111001;
				seg4 = 7'b0110000;
			end 
			32: begin 
				seg3 = 7'b0100100;
				seg4 = 7'b0110000;
			end 
			33: begin 
				seg3 = 7'b0110000;
				seg4 = 7'b0110000;
			end 
			34: begin 
				seg3 = 7'b0011001;
				seg4 = 7'b0110000;
			end 
			35: begin
				seg3 = 7'b0010010;
				seg4 = 7'b0110000;
			end 
			36: begin 
				seg3 = 7'b0000010;
				seg4 = 7'b0110000;
			end
			37: begin 
				seg3 = 7'b1111000;
				seg4 = 7'b0110000;
			end 
			38: begin 
				seg3 = 7'b0000000;
				seg4 = 7'b0110000;
			end
			39: begin 
				seg3 = 7'b0010000;
				seg4 = 7'b0110000;
			end 
			40: begin 
				seg3 = 7'b1000000;
				seg4 = 7'b0011001;
			end
			41: begin 
				seg3 = 7'b1111001;
				seg4 = 7'b0011001;
			end
			42: begin 
				seg3 = 7'b0100100;
				seg4 = 7'b0011001;
			end
			43: begin 
				seg3 = 7'b0110000;
				seg4 = 7'b0011001;
			end
			44: begin 
				seg3 = 7'b0011001;
				seg4 = 7'b0011001;
			end 
			45: begin 
				seg3 = 7'b0010010;
				seg4 = 7'b0011001;
			end 
			46: begin 
				seg3 = 7'b0000010;
				seg4 = 7'b0011001;
			end
			47: begin 
				seg3 = 7'b1111000;
				seg4 = 7'b0011001;
			end 
			48: begin 
				seg3 = 7'b0000000;
				seg4 = 7'b0011001;
			end 
			49: begin 
				seg3 = 7'b0010000;
				seg4 = 7'b0011001;
			end 
			50: begin 
				seg3 = 7'b1000000;
				seg4 = 7'b0010010;
			end
			51: begin 
				seg3 = 7'b1111001;
				seg4 = 7'b0010010;
			end
			52: begin 
				seg3 = 7'b0100100;
				seg4 = 7'b0010010;
			end
			53: begin 
				seg3 = 7'b0110000;
				seg4 = 7'b0010010;
			end
			54: begin 
				seg3 = 7'b0011001;
				seg4 = 7'b0010010;
			end 
			55: begin 
				seg3 = 7'b0010010;
				seg4 = 7'b0010010;
			end 
			56: begin 
				seg3 = 7'b0000010;
				seg4 = 7'b0010010;
			end
			57: begin 
				seg3 = 7'b1111000;
				seg4 = 7'b0010010;
			end 
			58: begin 
				seg3 = 7'b0000000;
				seg4 = 7'b0010010;
			end 
			59: begin 
				seg3 = 7'b0010000;
				seg4 = 7'b0010010;
			end 
			60: begin 
				seg3 = 7'b1000000;
				seg4 = 7'b0000010;
			end
			61: begin 
				seg3 = 7'b1111001;
				seg4 = 7'b0000010;
			end
			62: begin 
				seg3 = 7'b0100100;
				seg4 = 7'b0000010;
			end
			63: begin 
				seg3 = 7'b0110000;
				seg4 = 7'b0000010;
			end
			64: begin 
				seg3 = 7'b0011001;
				seg4 = 7'b0000010;
			end
		endcase
	
	end 

endmodule

