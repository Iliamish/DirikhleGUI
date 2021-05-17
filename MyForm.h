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
	using namespace ZedGraph;

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
	private: ZedGraph::ZedGraphControl^ zedGraphControl1;
	private: System::Windows::Forms::TabControl^ tabControl1;
	private: System::Windows::Forms::TabPage^ tabPage1;
	private: System::Windows::Forms::TabPage^ tabPage2;
	private: ZedGraph::ZedGraphControl^ zedGraphControl2;
	private: System::Windows::Forms::TabPage^ tabPage3;
	private: System::Windows::Forms::RadioButton^ radioButton1;
	private: System::Windows::Forms::RadioButton^ radioButton2;
	private: System::Windows::Forms::RadioButton^ radioButton3;
	private: ZedGraph::ZedGraphControl^ zedGraphControl3;
	private: System::Windows::Forms::TabControl^ tabControl2;
	private: System::Windows::Forms::TabPage^ tabPage4;
	private: System::Windows::Forms::TabPage^ tabPage5;
	private: System::Windows::Forms::DataGridView^ dataGridView2;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ dataGridViewTextBoxColumn1;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ dataGridViewTextBoxColumn2;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ dataGridViewTextBoxColumn3;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ dataGridViewTextBoxColumn4;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ dataGridViewTextBoxColumn5;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ dataGridViewTextBoxColumn6;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ dataGridViewTextBoxColumn7;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ Column3;
	private: System::Windows::Forms::TabPage^ tabPage6;
	private: System::Windows::Forms::TabPage^ tabPage7;
	private: System::Windows::Forms::TabPage^ tabPage8;
	private: ZedGraph::ZedGraphControl^ zedGraphControl4;
	private: ZedGraph::ZedGraphControl^ zedGraphControl5;
	private: ZedGraph::ZedGraphControl^ zedGraphControl6;












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
			this->components = (gcnew System::ComponentModel::Container());
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
			this->zedGraphControl1 = (gcnew ZedGraph::ZedGraphControl());
			this->tabControl1 = (gcnew System::Windows::Forms::TabControl());
			this->tabPage1 = (gcnew System::Windows::Forms::TabPage());
			this->tabPage2 = (gcnew System::Windows::Forms::TabPage());
			this->zedGraphControl2 = (gcnew ZedGraph::ZedGraphControl());
			this->tabPage3 = (gcnew System::Windows::Forms::TabPage());
			this->zedGraphControl3 = (gcnew ZedGraph::ZedGraphControl());
			this->tabPage6 = (gcnew System::Windows::Forms::TabPage());
			this->tabPage7 = (gcnew System::Windows::Forms::TabPage());
			this->tabPage8 = (gcnew System::Windows::Forms::TabPage());
			this->radioButton1 = (gcnew System::Windows::Forms::RadioButton());
			this->radioButton2 = (gcnew System::Windows::Forms::RadioButton());
			this->radioButton3 = (gcnew System::Windows::Forms::RadioButton());
			this->tabControl2 = (gcnew System::Windows::Forms::TabControl());
			this->tabPage4 = (gcnew System::Windows::Forms::TabPage());
			this->tabPage5 = (gcnew System::Windows::Forms::TabPage());
			this->dataGridView2 = (gcnew System::Windows::Forms::DataGridView());
			this->dataGridViewTextBoxColumn1 = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->dataGridViewTextBoxColumn2 = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->dataGridViewTextBoxColumn3 = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->dataGridViewTextBoxColumn4 = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->dataGridViewTextBoxColumn5 = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->dataGridViewTextBoxColumn6 = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->dataGridViewTextBoxColumn7 = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->Column3 = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->zedGraphControl4 = (gcnew ZedGraph::ZedGraphControl());
			this->zedGraphControl5 = (gcnew ZedGraph::ZedGraphControl());
			this->zedGraphControl6 = (gcnew ZedGraph::ZedGraphControl());
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->dataGridView1))->BeginInit();
			this->tabControl1->SuspendLayout();
			this->tabPage1->SuspendLayout();
			this->tabPage2->SuspendLayout();
			this->tabPage3->SuspendLayout();
			this->tabPage6->SuspendLayout();
			this->tabPage7->SuspendLayout();
			this->tabPage8->SuspendLayout();
			this->tabControl2->SuspendLayout();
			this->tabPage4->SuspendLayout();
			this->tabPage5->SuspendLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->dataGridView2))->BeginInit();
			this->SuspendLayout();
			// 
			// label1
			// 
			this->label1->AutoSize = true;
			this->label1->Location = System::Drawing::Point(108, 52);
			this->label1->Name = L"label1";
			this->label1->Size = System::Drawing::Size(99, 13);
			this->label1->TabIndex = 0;
			this->label1->Text = L"Число разбиений ";
			this->label1->Click += gcnew System::EventHandler(this, &MyForm::label1_Click);
			// 
			// textBox1
			// 
			this->textBox1->Location = System::Drawing::Point(213, 52);
			this->textBox1->Name = L"textBox1";
			this->textBox1->Size = System::Drawing::Size(100, 20);
			this->textBox1->TabIndex = 1;
			this->textBox1->Text = L"10";
			// 
			// button1
			// 
			this->button1->Location = System::Drawing::Point(739, 555);
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
			this->dataGridView1->Location = System::Drawing::Point(-4, 0);
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
			this->label2->Location = System::Drawing::Point(22, 80);
			this->label2->Name = L"label2";
			this->label2->Size = System::Drawing::Size(185, 13);
			this->label2->TabIndex = 12;
			this->label2->Text = L"Кратность вспомогательной сетки";
			this->label2->Click += gcnew System::EventHandler(this, &MyForm::label2_Click);
			// 
			// textBox2
			// 
			this->textBox2->Location = System::Drawing::Point(213, 80);
			this->textBox2->Name = L"textBox2";
			this->textBox2->Size = System::Drawing::Size(100, 20);
			this->textBox2->TabIndex = 13;
			this->textBox2->Text = L"4";
			// 
			// zedGraphControl1
			// 
			this->zedGraphControl1->Location = System::Drawing::Point(-2, 0);
			this->zedGraphControl1->Name = L"zedGraphControl1";
			this->zedGraphControl1->ScrollGrace = 0;
			this->zedGraphControl1->ScrollMaxX = 0;
			this->zedGraphControl1->ScrollMaxY = 0;
			this->zedGraphControl1->ScrollMaxY2 = 0;
			this->zedGraphControl1->ScrollMinX = 0;
			this->zedGraphControl1->ScrollMinY = 0;
			this->zedGraphControl1->ScrollMinY2 = 0;
			this->zedGraphControl1->Size = System::Drawing::Size(631, 329);
			this->zedGraphControl1->TabIndex = 14;
			// 
			// tabControl1
			// 
			this->tabControl1->Controls->Add(this->tabPage1);
			this->tabControl1->Controls->Add(this->tabPage2);
			this->tabControl1->Controls->Add(this->tabPage3);
			this->tabControl1->Controls->Add(this->tabPage6);
			this->tabControl1->Controls->Add(this->tabPage7);
			this->tabControl1->Controls->Add(this->tabPage8);
			this->tabControl1->Location = System::Drawing::Point(21, 155);
			this->tabControl1->Name = L"tabControl1";
			this->tabControl1->SelectedIndex = 0;
			this->tabControl1->Size = System::Drawing::Size(637, 355);
			this->tabControl1->TabIndex = 15;
			// 
			// tabPage1
			// 
			this->tabPage1->Controls->Add(this->zedGraphControl1);
			this->tabPage1->Location = System::Drawing::Point(4, 22);
			this->tabPage1->Name = L"tabPage1";
			this->tabPage1->Padding = System::Windows::Forms::Padding(3);
			this->tabPage1->Size = System::Drawing::Size(629, 329);
			this->tabPage1->TabIndex = 0;
			this->tabPage1->Text = L"Функции";
			this->tabPage1->UseVisualStyleBackColor = true;
			// 
			// tabPage2
			// 
			this->tabPage2->Controls->Add(this->zedGraphControl2);
			this->tabPage2->Location = System::Drawing::Point(4, 22);
			this->tabPage2->Name = L"tabPage2";
			this->tabPage2->Padding = System::Windows::Forms::Padding(3);
			this->tabPage2->Size = System::Drawing::Size(629, 329);
			this->tabPage2->TabIndex = 1;
			this->tabPage2->Text = L"Первые производные";
			this->tabPage2->UseVisualStyleBackColor = true;
			// 
			// zedGraphControl2
			// 
			this->zedGraphControl2->Location = System::Drawing::Point(-4, 0);
			this->zedGraphControl2->Name = L"zedGraphControl2";
			this->zedGraphControl2->ScrollGrace = 0;
			this->zedGraphControl2->ScrollMaxX = 0;
			this->zedGraphControl2->ScrollMaxY = 0;
			this->zedGraphControl2->ScrollMaxY2 = 0;
			this->zedGraphControl2->ScrollMinX = 0;
			this->zedGraphControl2->ScrollMinY = 0;
			this->zedGraphControl2->ScrollMinY2 = 0;
			this->zedGraphControl2->Size = System::Drawing::Size(631, 329);
			this->zedGraphControl2->TabIndex = 15;
			this->zedGraphControl2->Load += gcnew System::EventHandler(this, &MyForm::zedGraphControl2_Load);
			// 
			// tabPage3
			// 
			this->tabPage3->Controls->Add(this->zedGraphControl3);
			this->tabPage3->Location = System::Drawing::Point(4, 22);
			this->tabPage3->Name = L"tabPage3";
			this->tabPage3->Padding = System::Windows::Forms::Padding(3);
			this->tabPage3->Size = System::Drawing::Size(629, 329);
			this->tabPage3->TabIndex = 2;
			this->tabPage3->Text = L"Вторые производные";
			this->tabPage3->UseVisualStyleBackColor = true;
			// 
			// zedGraphControl3
			// 
			this->zedGraphControl3->Location = System::Drawing::Point(-4, 0);
			this->zedGraphControl3->Name = L"zedGraphControl3";
			this->zedGraphControl3->ScrollGrace = 0;
			this->zedGraphControl3->ScrollMaxX = 0;
			this->zedGraphControl3->ScrollMaxY = 0;
			this->zedGraphControl3->ScrollMaxY2 = 0;
			this->zedGraphControl3->ScrollMinX = 0;
			this->zedGraphControl3->ScrollMinY = 0;
			this->zedGraphControl3->ScrollMinY2 = 0;
			this->zedGraphControl3->Size = System::Drawing::Size(631, 329);
			this->zedGraphControl3->TabIndex = 16;
			// 
			// tabPage6
			// 
			this->tabPage6->Controls->Add(this->zedGraphControl4);
			this->tabPage6->Location = System::Drawing::Point(4, 22);
			this->tabPage6->Name = L"tabPage6";
			this->tabPage6->Size = System::Drawing::Size(629, 329);
			this->tabPage6->TabIndex = 3;
			this->tabPage6->Text = L"Разность функций";
			this->tabPage6->UseVisualStyleBackColor = true;
			// 
			// tabPage7
			// 
			this->tabPage7->Controls->Add(this->zedGraphControl5);
			this->tabPage7->Location = System::Drawing::Point(4, 22);
			this->tabPage7->Name = L"tabPage7";
			this->tabPage7->Padding = System::Windows::Forms::Padding(3);
			this->tabPage7->Size = System::Drawing::Size(629, 329);
			this->tabPage7->TabIndex = 4;
			this->tabPage7->Text = L"Разность первых производных";
			this->tabPage7->UseVisualStyleBackColor = true;
			// 
			// tabPage8
			// 
			this->tabPage8->Controls->Add(this->zedGraphControl6);
			this->tabPage8->Location = System::Drawing::Point(4, 22);
			this->tabPage8->Name = L"tabPage8";
			this->tabPage8->Padding = System::Windows::Forms::Padding(3);
			this->tabPage8->Size = System::Drawing::Size(629, 329);
			this->tabPage8->TabIndex = 5;
			this->tabPage8->Text = L"Разность вторых производных";
			this->tabPage8->UseVisualStyleBackColor = true;
			// 
			// radioButton1
			// 
			this->radioButton1->AutoSize = true;
			this->radioButton1->Checked = true;
			this->radioButton1->Location = System::Drawing::Point(352, 37);
			this->radioButton1->Name = L"radioButton1";
			this->radioButton1->Size = System::Drawing::Size(111, 17);
			this->radioButton1->TabIndex = 16;
			this->radioButton1->TabStop = true;
			this->radioButton1->Text = L"Тестовая задача";
			this->radioButton1->UseVisualStyleBackColor = true;
			this->radioButton1->CheckedChanged += gcnew System::EventHandler(this, &MyForm::radioButton1_CheckedChanged);
			// 
			// radioButton2
			// 
			this->radioButton2->AutoSize = true;
			this->radioButton2->Location = System::Drawing::Point(352, 94);
			this->radioButton2->Name = L"radioButton2";
			this->radioButton2->Size = System::Drawing::Size(122, 17);
			this->radioButton2->TabIndex = 17;
			this->radioButton2->Text = L"Основная задача 2";
			this->radioButton2->UseVisualStyleBackColor = true;
			// 
			// radioButton3
			// 
			this->radioButton3->AutoSize = true;
			this->radioButton3->Location = System::Drawing::Point(352, 65);
			this->radioButton3->Name = L"radioButton3";
			this->radioButton3->Size = System::Drawing::Size(122, 17);
			this->radioButton3->TabIndex = 18;
			this->radioButton3->Text = L"Основная задача 1";
			this->radioButton3->UseVisualStyleBackColor = true;
			this->radioButton3->CheckedChanged += gcnew System::EventHandler(this, &MyForm::radioButton3_CheckedChanged);
			// 
			// tabControl2
			// 
			this->tabControl2->Controls->Add(this->tabPage4);
			this->tabControl2->Controls->Add(this->tabPage5);
			this->tabControl2->Location = System::Drawing::Point(739, 155);
			this->tabControl2->Name = L"tabControl2";
			this->tabControl2->SelectedIndex = 0;
			this->tabControl2->Size = System::Drawing::Size(677, 355);
			this->tabControl2->TabIndex = 19;
			// 
			// tabPage4
			// 
			this->tabPage4->Controls->Add(this->dataGridView1);
			this->tabPage4->Location = System::Drawing::Point(4, 22);
			this->tabPage4->Name = L"tabPage4";
			this->tabPage4->Padding = System::Windows::Forms::Padding(3);
			this->tabPage4->Size = System::Drawing::Size(669, 329);
			this->tabPage4->TabIndex = 0;
			this->tabPage4->Text = L"Коэффициенты сплайна";
			this->tabPage4->UseVisualStyleBackColor = true;
			// 
			// tabPage5
			// 
			this->tabPage5->Controls->Add(this->dataGridView2);
			this->tabPage5->Location = System::Drawing::Point(4, 22);
			this->tabPage5->Name = L"tabPage5";
			this->tabPage5->Padding = System::Windows::Forms::Padding(3);
			this->tabPage5->Size = System::Drawing::Size(669, 329);
			this->tabPage5->TabIndex = 1;
			this->tabPage5->Text = L"Результаты";
			this->tabPage5->UseVisualStyleBackColor = true;
			// 
			// dataGridView2
			// 
			this->dataGridView2->ColumnHeadersHeightSizeMode = System::Windows::Forms::DataGridViewColumnHeadersHeightSizeMode::AutoSize;
			this->dataGridView2->Columns->AddRange(gcnew cli::array< System::Windows::Forms::DataGridViewColumn^  >(8) {
				this->dataGridViewTextBoxColumn1,
					this->dataGridViewTextBoxColumn2, this->dataGridViewTextBoxColumn3, this->dataGridViewTextBoxColumn4, this->dataGridViewTextBoxColumn5,
					this->dataGridViewTextBoxColumn6, this->dataGridViewTextBoxColumn7, this->Column3
			});
			this->dataGridView2->Location = System::Drawing::Point(-4, 0);
			this->dataGridView2->Name = L"dataGridView2";
			this->dataGridView2->RowHeadersVisible = false;
			this->dataGridView2->Size = System::Drawing::Size(677, 333);
			this->dataGridView2->TabIndex = 12;
			// 
			// dataGridViewTextBoxColumn1
			// 
			this->dataGridViewTextBoxColumn1->HeaderText = L"j";
			this->dataGridViewTextBoxColumn1->Name = L"dataGridViewTextBoxColumn1";
			this->dataGridViewTextBoxColumn1->ReadOnly = true;
			this->dataGridViewTextBoxColumn1->Width = 50;
			// 
			// dataGridViewTextBoxColumn2
			// 
			this->dataGridViewTextBoxColumn2->HeaderText = L"Xj";
			this->dataGridViewTextBoxColumn2->Name = L"dataGridViewTextBoxColumn2";
			this->dataGridViewTextBoxColumn2->ReadOnly = true;
			this->dataGridViewTextBoxColumn2->Width = 50;
			// 
			// dataGridViewTextBoxColumn3
			// 
			this->dataGridViewTextBoxColumn3->HeaderText = L"F (xj)";
			this->dataGridViewTextBoxColumn3->Name = L"dataGridViewTextBoxColumn3";
			this->dataGridViewTextBoxColumn3->ReadOnly = true;
			this->dataGridViewTextBoxColumn3->Width = 90;
			// 
			// dataGridViewTextBoxColumn4
			// 
			this->dataGridViewTextBoxColumn4->HeaderText = L"S (xj)";
			this->dataGridViewTextBoxColumn4->Name = L"dataGridViewTextBoxColumn4";
			this->dataGridViewTextBoxColumn4->ReadOnly = true;
			this->dataGridViewTextBoxColumn4->Width = 90;
			// 
			// dataGridViewTextBoxColumn5
			// 
			this->dataGridViewTextBoxColumn5->HeaderText = L"F (xj) - S(xj)";
			this->dataGridViewTextBoxColumn5->Name = L"dataGridViewTextBoxColumn5";
			this->dataGridViewTextBoxColumn5->ReadOnly = true;
			// 
			// dataGridViewTextBoxColumn6
			// 
			this->dataGridViewTextBoxColumn6->HeaderText = L"F \' (xj)";
			this->dataGridViewTextBoxColumn6->Name = L"dataGridViewTextBoxColumn6";
			// 
			// dataGridViewTextBoxColumn7
			// 
			this->dataGridViewTextBoxColumn7->HeaderText = L"S\' (xj)";
			this->dataGridViewTextBoxColumn7->Name = L"dataGridViewTextBoxColumn7";
			// 
			// Column3
			// 
			this->Column3->HeaderText = L"F \' (xj) - S\'(xj)";
			this->Column3->Name = L"Column3";
			// 
			// zedGraphControl4
			// 
			this->zedGraphControl4->Location = System::Drawing::Point(-2, 0);
			this->zedGraphControl4->Name = L"zedGraphControl4";
			this->zedGraphControl4->ScrollGrace = 0;
			this->zedGraphControl4->ScrollMaxX = 0;
			this->zedGraphControl4->ScrollMaxY = 0;
			this->zedGraphControl4->ScrollMaxY2 = 0;
			this->zedGraphControl4->ScrollMinX = 0;
			this->zedGraphControl4->ScrollMinY = 0;
			this->zedGraphControl4->ScrollMinY2 = 0;
			this->zedGraphControl4->Size = System::Drawing::Size(631, 329);
			this->zedGraphControl4->TabIndex = 17;
			// 
			// zedGraphControl5
			// 
			this->zedGraphControl5->Location = System::Drawing::Point(-4, 0);
			this->zedGraphControl5->Name = L"zedGraphControl5";
			this->zedGraphControl5->ScrollGrace = 0;
			this->zedGraphControl5->ScrollMaxX = 0;
			this->zedGraphControl5->ScrollMaxY = 0;
			this->zedGraphControl5->ScrollMaxY2 = 0;
			this->zedGraphControl5->ScrollMinX = 0;
			this->zedGraphControl5->ScrollMinY = 0;
			this->zedGraphControl5->ScrollMinY2 = 0;
			this->zedGraphControl5->Size = System::Drawing::Size(631, 329);
			this->zedGraphControl5->TabIndex = 18;
			// 
			// zedGraphControl6
			// 
			this->zedGraphControl6->Location = System::Drawing::Point(-4, 0);
			this->zedGraphControl6->Name = L"zedGraphControl6";
			this->zedGraphControl6->ScrollGrace = 0;
			this->zedGraphControl6->ScrollMaxX = 0;
			this->zedGraphControl6->ScrollMaxY = 0;
			this->zedGraphControl6->ScrollMaxY2 = 0;
			this->zedGraphControl6->ScrollMinX = 0;
			this->zedGraphControl6->ScrollMinY = 0;
			this->zedGraphControl6->ScrollMinY2 = 0;
			this->zedGraphControl6->Size = System::Drawing::Size(631, 329);
			this->zedGraphControl6->TabIndex = 19;
			// 
			// MyForm
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(6, 13);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->ClientSize = System::Drawing::Size(1540, 648);
			this->Controls->Add(this->tabControl2);
			this->Controls->Add(this->radioButton3);
			this->Controls->Add(this->radioButton2);
			this->Controls->Add(this->radioButton1);
			this->Controls->Add(this->tabControl1);
			this->Controls->Add(this->textBox2);
			this->Controls->Add(this->label2);
			this->Controls->Add(this->label7);
			this->Controls->Add(this->button1);
			this->Controls->Add(this->textBox1);
			this->Controls->Add(this->label1);
			this->Name = L"MyForm";
			this->Text = L"MyForm";
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->dataGridView1))->EndInit();
			this->tabControl1->ResumeLayout(false);
			this->tabPage1->ResumeLayout(false);
			this->tabPage2->ResumeLayout(false);
			this->tabPage3->ResumeLayout(false);
			this->tabPage6->ResumeLayout(false);
			this->tabPage7->ResumeLayout(false);
			this->tabPage8->ResumeLayout(false);
			this->tabControl2->ResumeLayout(false);
			this->tabPage4->ResumeLayout(false);
			this->tabPage5->ResumeLayout(false);
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->dataGridView2))->EndInit();
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

		GraphPane^ panel = zedGraphControl1->GraphPane;
		GraphPane^ panel2 = zedGraphControl2->GraphPane;
		GraphPane^ panel3 = zedGraphControl3->GraphPane;
		GraphPane^ panel4 = zedGraphControl4->GraphPane;
		GraphPane^ panel5 = zedGraphControl5->GraphPane;
		GraphPane^ panel6 = zedGraphControl6->GraphPane;

		panel3->CurveList->Clear();
		panel->CurveList->Clear();
		panel2->CurveList->Clear();
		double fa = 0;
		double fb = 0;
		int a, b;
		int n= Convert::ToDouble(textBox1->Text); // максимальное число итераций (не менее 1)
		int N = Convert::ToDouble(textBox2->Text);
		std::function<double(double)> fun;
		std::function<double(double)> fun1;
		std::function<double(double)> fun2;
	
		if (radioButton1->Checked) {
		   fun = testFunc;
		    fun1 = testFuncDiv1;
			fun2 = testFuncDiv2;
			a = -1;
			b = 1;
			panel->XAxis->Scale->Min = -1.5;
			panel->XAxis->Scale->Max = 1.5;

			panel2->XAxis->Scale->Min = -1.5;
			panel2->XAxis->Scale->Max = 1.5;

			panel3->XAxis->Scale->Min = -1.5;
			panel3->XAxis->Scale->Max = 1.5;


		}
		else if (radioButton3->Checked) {
			 fun =mainFunc;
			 fun1 = mainFuncDiv1;
			 fun2 = mainFuncDiv2;
			a = 0;
			b = 1;
			panel->XAxis->Scale->Min = -0.1;
			panel->XAxis->Scale->Max = 1.1;

			panel2->XAxis->Scale->Min = -0.1;
			panel2->XAxis->Scale->Max = 1.1;

			panel3->XAxis->Scale->Min = -0.1;
			panel3->XAxis->Scale->Max = 1.1;
		}
		else {
			a = 0;
			b = 1;
			fun = mainFunc1;
			 fun1 = mainFunc1Div1;
			 fun2 = mainFunc1Div2;
		}
		double f_x_max = 0;
		double f_max = 0;

		double f1_x_max = 0;
		double f1_max = 0;

		double f2_x_max = 0;
		double f2_max = 0;

		
		PointPairList^ f_list = gcnew ZedGraph::PointPairList();
		PointPairList^ f1_list = gcnew ZedGraph::PointPairList();
		PointPairList^ f2_list = gcnew ZedGraph::PointPairList();

		PointPairList^ s_list = gcnew ZedGraph::PointPairList();
		PointPairList^ s1_list = gcnew ZedGraph::PointPairList();
		PointPairList^ s2_list = gcnew ZedGraph::PointPairList();

		PointPairList^ dif = gcnew ZedGraph::PointPairList();
		PointPairList^ dif1 = gcnew ZedGraph::PointPairList();
		PointPairList^ dif2 = gcnew ZedGraph::PointPairList();

		TResults result= tfunc1(a, b, n, fun);
		double h = double(b - a) / n;
		std::vector<double> x;
		std::vector<double> f;
		std::vector<double> f1;
		std::vector<double> f2;
		
		//fa = mainFunc(a);
	
		

		double f_curr = 0;
		double current_x = a;
		double h_n = double(b-a) / (N * n);
		f.push_back(result.a[0] + result.b[0] * (current_x - result.res_vec[0 + 1].first) + (result.c[0 + 1] / 2) * pow(current_x - result.res_vec[0 + 1].first, 2) + (result.d[0] / 6) * pow(current_x - result.res_vec[0 + 1].first, 3));
		f1.push_back(result.b[0] + (result.c[0 + 1]) * pow(current_x - result.res_vec[0 + 1].first, 1) + (result.d[0] / 2) * pow(current_x - result.res_vec[0 + 1].first, 2));
		f2.push_back(result.c[0 + 1] + (result.d[0]) * pow(current_x - result.res_vec[0 + 1].first, 1));
		x.push_back(a);
		dif->Add(current_x, abs(fun(current_x) - f.back()));
		dif1->Add(current_x, abs(fun1(current_x) - f1.back()));
		dif2->Add(current_x, abs(fun2(current_x) - f2.back()));
		int k = 0;
		dataGridView2->Rows->Add();
		dataGridView2->Rows[0]->Cells[0]->Value = k;
		dataGridView2->Rows[0]->Cells[1]->Value =a;
		dataGridView2->Rows[0]->Cells[2]->Value =fun(a) ;
		dataGridView2->Rows[0]->Cells[3]->Value = f.back();
		dataGridView2->Rows[0]->Cells[4]->Value = abs(fun(a) - f.back());
		dataGridView2->Rows[0]->Cells[5]->Value = fun1(a);
		dataGridView2->Rows[0]->Cells[6]->Value = f1.back();
		dataGridView2->Rows[0]->Cells[7]->Value = abs(fun1(current_x) - f1.back());

		for (size_t i = 0; i < result.b.size(); i++) {
			//x.push_back(result.res_vec[i].first);
			for (int j = 1; j < N; j++) {
				current_x += h_n;
				x.push_back(current_x);
				k++;
				// Подсчет для основной функции
				f_curr = result.a[i] + result.b[i] * (current_x - result.res_vec[i + 1].first) + (result.c[i + 1] / 2) * pow(current_x - result.res_vec[i + 1].first, 2) + (result.d[i] / 6) * pow(current_x - result.res_vec[i + 1].first, 3);
				f.push_back(f_curr);
				if (abs(fun(current_x) - f_curr)>f_max) {
					f_max = abs(fun(current_x) - f_curr);
					f_x_max = current_x;
				}

				// Подсчет для первой производной
				f_curr = result.b[i] + (result.c[i + 1]) * pow(current_x - result.res_vec[i + 1].first, 1) + (result.d[i] / 2) * pow(current_x - result.res_vec[i + 1].first, 2);
				f1.push_back(f_curr);
				if (abs(fun1(current_x)-f_curr)>f1_max) {
					f1_max = abs(fun1(current_x) - f_curr);
					f_x_max = current_x;
				}

				// Подсчет для второй производной
				f_curr = result.c[i+1]+ (result.d[i] ) * pow(current_x - result.res_vec[i + 1].first, 1);
				if (abs(fun2(current_x) - f_curr) > f2_max) {
					f2_max = abs(fun2(current_x) - f_curr);
					f_x_max = current_x;
				}
				dif->Add(current_x, abs(fun(current_x) - f.back()));
				dif1->Add(current_x, abs(fun1(current_x) - f1.back()));
				dif2->Add(current_x, abs(fun2(current_x) - f2.back()));

				double res_func = fun(current_x);
				double res_func1 = fun1(current_x);
				dataGridView2->Rows->Add();
				dataGridView2->Rows[k]->Cells[0]->Value = k;
				dataGridView2->Rows[k]->Cells[1]->Value = current_x;
				dataGridView2->Rows[k]->Cells[2]->Value = res_func;
				dataGridView2->Rows[k]->Cells[3]->Value = f.back();
				dataGridView2->Rows[k]->Cells[4]->Value = abs(res_func - f.back());
				dataGridView2->Rows[k]->Cells[5]->Value = res_func1;
				dataGridView2->Rows[k]->Cells[6]->Value = f1.back();
				dataGridView2->Rows[k]->Cells[7]->Value = abs(res_func1- f1.back());
			}
			k++;
			f.push_back(result.a[i]);
			f1.push_back(result.b[i]);
			f2.push_back(result.c[i + 1]);
			current_x += h_n;

			dif->Add(current_x, abs(fun(current_x) - f.back()));
			dif1->Add(current_x, abs(fun1(current_x) - f1.back()));
			dif2->Add(current_x, abs(fun2(current_x) - f2.back()));

			double res_func = fun(current_x);
			double res_func1 = fun1(current_x);
			dataGridView2->Rows->Add();
			dataGridView2->Rows[k]->Cells[0]->Value = k;
			dataGridView2->Rows[k]->Cells[1]->Value = current_x;
			dataGridView2->Rows[k]->Cells[2]->Value =res_func;
			dataGridView2->Rows[k]->Cells[3]->Value = f.back();
			dataGridView2->Rows[k]->Cells[4]->Value = equalsZero(abs(res_func - f.back()));
			dataGridView2->Rows[k]->Cells[5]->Value =res_func1;
			dataGridView2->Rows[k]->Cells[6]->Value = f1.back();
			dataGridView2->Rows[k]->Cells[7]->Value = abs(res_func1 - f1.back());
			if (abs(fun(result.res_vec[i + 1].first) -result.a[i]) > f_max) {
				f_max = abs(fun(result.res_vec[i + 1].first) - result.a[i]);
				f_x_max = result.res_vec[i + 1].first;
			}
			if (abs(fun1(result.res_vec[i + 1].first) - result.b[i]) > f1_max) {
				f1_max = abs(fun1(result.res_vec[i + 1].first) - result.b[i]);
				f1_x_max = result.res_vec[i + 1].first;
			}
			if (abs(fun2(result.res_vec[i + 1].first) - result.c[i+1]) > f2_max) {
				f2_max = abs(fun2(result.res_vec[i + 1].first) - result.c[i+1]);
				f2_x_max = result.res_vec[i + 1].first;
			}
			x.push_back(result.res_vec[i + 1].first);

			dataGridView1->Rows->Add();
			dataGridView1->Rows[i]->Cells[0]->Value = i+1;
			dataGridView1->Rows[i]->Cells[1]->Value =round(1000*result.res_vec[i].first)/1000;
			dataGridView1->Rows[i]->Cells[2]->Value =round(1000*result.res_vec[i + 1].first)/1000;
			dataGridView1->Rows[i]->Cells[3]->Value = round(1000*result.a[i])/1000;
			dataGridView1->Rows[i]->Cells[4]->Value = round(1000*result.b[i])/1000;
			dataGridView1->Rows[i]->Cells[5]->Value = round(1000*result.c[i+1])/1000;
			dataGridView1->Rows[i]->Cells[6]->Value = round (1000*result.d[i])/1000;
		}
		double step = double(b - a) / 100;
		current_x = a;
		int i = 0;
		s_list->Add(a, fun(a));
		f_list->Add(a, fun(a));

		s1_list->Add(a, result.b[0] + (result.c[0 + 1]) * pow(current_x - result.res_vec[0 + 1].first, 1) + (result.d[0] / 2) * pow(current_x - result.res_vec[0 + 1].first, 2));
		f1_list->Add(a, fun1(a));

		s2_list->Add(a, result.c[0 + 1] + (result.d[0]) * pow(current_x - result.res_vec[0 + 1].first, 1));
		f2_list->Add(a, fun2(a));

		for (int j = 0; j < 99; j++) {
			current_x += step;
			while (current_x > result.res_vec[i + 1].first) {
				i++;
			}
			f_curr = result.a[i] + result.b[i] * (current_x - result.res_vec[i + 1].first) + (result.c[i + 1] / 2) * pow(current_x - result.res_vec[i + 1].first, 2) + (result.d[i] / 6) * pow(current_x - result.res_vec[i + 1].first, 3);
			s_list->Add(current_x, f_curr);
			f_list->Add(current_x, fun(current_x));

			f_curr = result.b[i] + (result.c[i + 1]) * pow(current_x - result.res_vec[i + 1].first, 1) + (result.d[i] / 2) * pow(current_x - result.res_vec[i + 1].first, 2);
			s1_list->Add(current_x, f_curr);
			f1_list->Add(current_x, fun1(current_x));

			f_curr = result.c[i + 1] + (result.d[i]) * pow(current_x - result.res_vec[i + 1].first, 1);
			s2_list->Add(current_x, f_curr);
			f2_list->Add(current_x, fun2(current_x));

		}
		s_list->Add(b, fun(b));
		f_list->Add(b, fun(b));

		s1_list->Add(b, result.b[i]);
		f1_list->Add(b, fun1(b));

		s2_list->Add(b, result.c[i + 1] );
		f2_list->Add(b, fun2(b));

		LineItem Curve1 = panel->AddCurve("S ( x )", s_list, Color::Red, SymbolType::None);
		LineItem Curve2 = panel->AddCurve("f ( x )", f_list, Color::Blue, SymbolType::None);

		LineItem Curve3 = panel2->AddCurve("S' ( x )", s1_list, Color::Red, SymbolType::None);
		LineItem Curve4 = panel2->AddCurve("f' ( x )", f1_list, Color::Blue, SymbolType::None);

		LineItem Curve5 = panel3->AddCurve("S'' ( x )", s2_list, Color::Red, SymbolType::None);
		LineItem Curve6 = panel3->AddCurve("f'' ( x )", f2_list, Color::Blue, SymbolType::None);
		
		LineItem Curve7 = panel4->AddCurve("|F(x) - S(x)|", dif, Color::Red, SymbolType::None);
		LineItem Curve8 = panel5->AddCurve("|F'(x) - S'(x)|", dif1, Color::Blue, SymbolType::None);
		LineItem Curve9 = panel6->AddCurve("|F''(x)- S''(x)|", dif2, Color::Blue, SymbolType::None);

		zedGraphControl1->AxisChange();
		// Обновляем график
		zedGraphControl1->Invalidate();

		zedGraphControl2->AxisChange();
		// Обновляем график
		zedGraphControl2->Invalidate();

		zedGraphControl3->AxisChange();
		// Обновляем график
		zedGraphControl3->Invalidate();
		String ^string1 =gcnew String("                                                                 Справка\n                                                      Сетка сплайна: n = "+n +"\n                                                  Контрольная сетка: N = "+N*n+"\n                               Погрешность сплайна на контрольной сетке"+"\n              max |F(x) - S(x)| = "+f_max +" при x = "+f_x_max+"\n                           Погрешность производной на контрольной сетке"+"\n                     max|F'(x) - S'(x)| = "+f1_max +" при x = "+ f1_x_max +"\n                Погрешность второй производной на контрольной сетке "+ "\n                             max |F''(x) - S''(x)| = "+f2_max+"при x = "+ f2_x_max);
		MessageBox::Show(string1);

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
private: System::Void zedGraphControl2_Load(System::Object^ sender, System::EventArgs^ e) {
}
private: System::Void radioButton1_CheckedChanged(System::Object^ sender, System::EventArgs^ e) {
}
private: System::Void radioButton3_CheckedChanged(System::Object^ sender, System::EventArgs^ e) {
}
};
}
