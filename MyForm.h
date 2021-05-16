#pragma once
#include "Dirikhle.h"
#include <thread>
#include <msclr\marshal_cppstd.h>

namespace DirikhleGUI {

	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;
	using namespace System::Threading;
	using namespace System::Diagnostics;
	using namespace System::ComponentModel;

	/// <summary>
	/// Summary for MyForm
	/// </summary>
	public ref class MyForm : public System::Windows::Forms::Form
	{
	public:
		MyForm(void)
		{
			InitializeComponent();
			//
			//TODO: Add the constructor code here
			//
		}

	protected:
		/// <summary>
		/// Clean up any resources being used.
		/// </summary>
		~MyForm()
		{
			if (components)
			{
				delete components;
			}
		}
	private: System::Windows::Forms::Label^ label1;
	protected:
	private: System::Windows::Forms::TextBox^ textBox1;










	private: System::Windows::Forms::Button^ button1;
	private: System::ComponentModel::IContainer^ components;

	private: System::Windows::Forms::Label^ label7;




	private: System::Windows::Forms::DataGridView^ dataGridView1;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ NumberOfColumn;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ X;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ V;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ UandDouble;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ SubVandSome;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ Column1;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ Column2;
	private: System::Windows::Forms::Label^ label2;
	private: System::Windows::Forms::TextBox^ textBox2;











	private: Thread^ thread2 = nullptr;

	private:
		/// <summary>
		/// Required designer variable.
		/// </summary>


#pragma region Windows Form Designer generated code
		/// <summary>
		/// Required method for Designer support - do not modify
		/// the contents of this method with the code editor.
		/// </summary>
		void InitializeComponent(void)
		{
			this->label1 = (gcnew System::Windows::Forms::Label());
			this->textBox1 = (gcnew System::Windows::Forms::TextBox());
			this->button1 = (gcnew System::Windows::Forms::Button());
			this->label7 = (gcnew System::Windows::Forms::Label());
			this->dataGridView1 = (gcnew System::Windows::Forms::DataGridView());
			this->NumberOfColumn = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->X = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->V = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->UandDouble = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->SubVandSome = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->Column1 = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->Column2 = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->label2 = (gcnew System::Windows::Forms::Label());
			this->textBox2 = (gcnew System::Windows::Forms::TextBox());
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->dataGridView1))->BeginInit();
			this->SuspendLayout();
			// 
			// label1
			// 
			this->label1->AutoSize = true;
			this->label1->Location = System::Drawing::Point(106, 86);
			this->label1->Name = L"label1";
			this->label1->Size = System::Drawing::Size(99, 13);
			this->label1->TabIndex = 0;
			this->label1->Text = L"Число разбиений ";
			this->label1->Click += gcnew System::EventHandler(this, &MyForm::label1_Click);
			// 
			// textBox1
			// 
			this->textBox1->Location = System::Drawing::Point(211, 86);
			this->textBox1->Name = L"textBox1";
			this->textBox1->Size = System::Drawing::Size(100, 20);
			this->textBox1->TabIndex = 1;
			this->textBox1->Text = L"10";
			// 
			// button1
			// 
			this->button1->Location = System::Drawing::Point(725, 12);
			this->button1->Name = L"button1";
			this->button1->Size = System::Drawing::Size(183, 46);
			this->button1->TabIndex = 4;
			this->button1->Text = L"Вычислить";
			this->button1->UseVisualStyleBackColor = true;
			this->button1->Click += gcnew System::EventHandler(this, &MyForm::button1_Click);
			// 
			// label7
			// 
			this->label7->AutoSize = true;
			this->label7->Location = System::Drawing::Point(15, 386);
			this->label7->Name = L"label7";
			this->label7->Size = System::Drawing::Size(0, 13);
			this->label7->TabIndex = 6;
			// 
			// dataGridView1
			// 
			this->dataGridView1->ColumnHeadersHeightSizeMode = System::Windows::Forms::DataGridViewColumnHeadersHeightSizeMode::AutoSize;
			this->dataGridView1->Columns->AddRange(gcnew cli::array< System::Windows::Forms::DataGridViewColumn^  >(7) {
				this->NumberOfColumn,
					this->X, this->V, this->UandDouble, this->SubVandSome, this->Column1, this->Column2
			});
			this->dataGridView1->Location = System::Drawing::Point(725, 86);
			this->dataGridView1->Name = L"dataGridView1";
			this->dataGridView1->RowHeadersVisible = false;
			this->dataGridView1->Size = System::Drawing::Size(677, 333);
			this->dataGridView1->TabIndex = 11;
			this->dataGridView1->CellContentClick += gcnew System::Windows::Forms::DataGridViewCellEventHandler(this, &MyForm::dataGridView1_CellContentClick);
			// 
			// NumberOfColumn
			// 
			this->NumberOfColumn->HeaderText = L"i";
			this->NumberOfColumn->Name = L"NumberOfColumn";
			this->NumberOfColumn->ReadOnly = true;
			// 
			// X
			// 
			this->X->HeaderText = L"Xi-1";
			this->X->Name = L"X";
			this->X->ReadOnly = true;
			this->X->Width = 50;
			// 
			// V
			// 
			this->V->HeaderText = L"Xi";
			this->V->Name = L"V";
			this->V->ReadOnly = true;
			// 
			// UandDouble
			// 
			this->UandDouble->HeaderText = L"ai";
			this->UandDouble->Name = L"UandDouble";
			this->UandDouble->ReadOnly = true;
			// 
			// SubVandSome
			// 
			this->SubVandSome->HeaderText = L"bi";
			this->SubVandSome->Name = L"SubVandSome";
			this->SubVandSome->ReadOnly = true;
			this->SubVandSome->Width = 124;
			// 
			// Column1
			// 
			this->Column1->HeaderText = L"ci";
			this->Column1->Name = L"Column1";
			// 
			// Column2
			// 
			this->Column2->HeaderText = L"di";
			this->Column2->Name = L"Column2";
			// 
			// label2
			// 
			this->label2->AutoSize = true;
			this->label2->Location = System::Drawing::Point(20, 114);
			this->label2->Name = L"label2";
			this->label2->Size = System::Drawing::Size(185, 13);
			this->label2->TabIndex = 12;
			this->label2->Text = L"Кратность вспомогательной сетки";
			this->label2->Click += gcnew System::EventHandler(this, &MyForm::label2_Click);
			// 
			// textBox2
			// 
			this->textBox2->Location = System::Drawing::Point(211, 114);
			this->textBox2->Name = L"textBox2";
			this->textBox2->Size = System::Drawing::Size(100, 20);
			this->textBox2->TabIndex = 13;
			this->textBox2->Text = L"4";
			// 
			// MyForm
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(6, 13);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->ClientSize = System::Drawing::Size(1540, 648);
			this->Controls->Add(this->textBox2);
			this->Controls->Add(this->label2);
			this->Controls->Add(this->dataGridView1);
			this->Controls->Add(this->label7);
			this->Controls->Add(this->button1);
			this->Controls->Add(this->textBox1);
			this->Controls->Add(this->label1);
			this->Name = L"MyForm";
			this->Text = L"MyForm";
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->dataGridView1))->EndInit();
			this->ResumeLayout(false);
			this->PerformLayout();

		}
#pragma endregion
		private: void BackgroundWorker1_DoWork() {
			// Sleep 2 seconds to emulate getting data.
			system("python show_plot.py");
		}


	private: System::Void pictureBox1_Click(System::Object^ sender, System::EventArgs^ e) {
	}


	private: System::Void button1_Click(System::Object^ sender, System::EventArgs^ e) {
		


		int n= Convert::ToDouble(textBox1->Text); // максимальное число итераций (не менее 1)
		int N = Convert::ToDouble(textBox2->Text);

		int a = -1;
		int b = 1;
		TResults result= tfunc1(a, b, n, testFunc);
	
		for (size_t i = 0; i < result.b.size(); i++) {

			dataGridView1->Rows->Add();
			dataGridView1->Rows[i]->Cells[0]->Value = i+1;
			dataGridView1->Rows[i]->Cells[1]->Value =round(1000*result.res_vec[i].first)/1000;
			dataGridView1->Rows[i]->Cells[2]->Value =round(1000*result.res_vec[i + 1].first)/1000;
			dataGridView1->Rows[i]->Cells[3]->Value = round(1000*result.a[i])/1000;
			dataGridView1->Rows[i]->Cells[4]->Value = round(1000*result.b[i])/1000;
			dataGridView1->Rows[i]->Cells[5]->Value = round(1000*result.c[i+1])/1000;
			dataGridView1->Rows[i]->Cells[6]->Value = round (1000*result.d[i])/1000;
		}
		
	}
private: System::Void checkBox1_CheckedChanged(System::Object^ sender, System::EventArgs^ e) {
}
private: System::Void listBox1_SelectedIndexChanged(System::Object^ sender, System::EventArgs^ e) {
}
private: System::Void radioButton2_CheckedChanged(System::Object^ sender, System::EventArgs^ e) {
}
private: System::Void dataGridView1_CellContentClick(System::Object^ sender, System::Windows::Forms::DataGridViewCellEventArgs^ e) {
}
private: System::Void label1_Click(System::Object^ sender, System::EventArgs^ e) {
}
private: System::Void label2_Click(System::Object^ sender, System::EventArgs^ e) {
}
};
}
